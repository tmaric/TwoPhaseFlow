#!/usr/bin/env python3
"""PML pulse verification for acousticWaveFoam.

Reads the physical-domain probes (p, pDot, U) and, if present, the
whole-domain probes (p, v = psi, u = w) and writes
postProcessing/pmlVerification/{summary.json,energy.csv,maxima.csv}:

- E2 of p and U against the free-space solution of a Gaussian pulse in a
  homogeneous gas (2D Hankel transform, 3D closed form), time-integrated
  over the samples up to T_obs; omitted when a drop is present;
- the physical energy E_ac (p, U) and the discrete wave energy E_h (q, p)
  of the physical domain, with energy weights V*beta and V*rho;
- the maxima of |p|, |psi| and |w| over the whole domain, PML included;
- growth rates over the last quarter of the run, classified as decaying,
  bounded or growing;
- the largest |p| in the physical domain after T_RESIDUAL, relative to P0.
  In 3D the exact pressure has vanished there once the pulse has passed,
  so this is PML reflection or a dispersive tail.
The comparison with the enlarged-domain reference is made per study by
summarize_pml_study.py, since it needs a second case.
"""
import json
import math
import warnings
from pathlib import Path

import numpy as np

CASE = Path(__file__).resolve().parent
OUT = CASE / "postProcessing" / "pmlVerification"
T_OBS = 4e-4
T_RESIDUAL = 3e-4
# A quantity whose growth factor over the last quarter stays within
# [1/BOUND, BOUND] is classified as bounded.
BOUND = 1.05


def read_probes(group, field, components=None):
    path = CASE / "postProcessing" / group / "0" / field
    times, rows = [], []
    with path.open() as stream:
        for line in stream:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.replace("(", " ").replace(")", " ").split()
            times.append(float(parts[0]))
            rows.append(np.array(parts[1:], dtype=float))
    # A run that blew up may end with a truncated row.
    if len(rows) > 1 and rows[-1].size != rows[0].size:
        times, rows = times[:-1], rows[:-1]
    data = np.array(rows)
    if components:
        data = data.reshape(len(times), -1, 3)[:, :, :components]
    return np.array(times), data


def free_space(s, points, times):
    """Exact pressure and velocity of the Gaussian pulse in free space."""
    p0, width, c = s["p0"], s["pulseWidth"], s["cG"]
    rho = s["rhoG"]
    d = points - np.array(s["source"])[None, :points.shape[1]]
    r = np.linalg.norm(d, axis=1)
    unit = d / np.where(r > 0, r, 1.0)[:, None]
    p = np.empty((len(times), len(r)))
    ur = np.empty_like(p)
    if s["dim"] == 2:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from scipy.special import j0, j1
        # Gauss-Legendre on k in [0, 9/s]; the weight exp(-k^2 s^2/2) is
        # below 1e-17 beyond. Tabulated on a fine radial grid and
        # interpolated to the cell centres.
        nodes, w = np.polynomial.legendre.leggauss(4000)
        k_max = 9.0 / width
        k = 0.5 * k_max * (nodes + 1.0)
        w = 0.5 * k_max * w * k * np.exp(-0.5 * (k * width) ** 2)
        r_grid = np.linspace(0.0, r.max() * (1 + 1e-9), 8001)
        kr = np.outer(r_grid, k)
        bessel0, bessel1 = j0(kr), j1(kr)
        for i, t in enumerate(times):
            p_line = p0 * width**2 * bessel0 @ (w * np.cos(c * k * t))
            u_line = (p0 * width**2 / (rho * c)
                      * bessel1 @ (w * np.sin(c * k * t)))
            p[i] = np.interp(r, r_grid, p_line)
            ur[i] = np.interp(r, r_grid, u_line)
    else:
        def g(xi):
            return np.exp(-0.5 * (xi / width) ** 2)
        for i, t in enumerate(times):
            ct = c * t
            rs = np.where(r > 1e-12, r, 1e-12)
            p[i] = p0 * ((rs - ct) * g(rs - ct) + (rs + ct) * g(rs + ct)) / (2 * rs)
        ur = None
    u = None if ur is None else ur[:, :, None] * unit[None, :, :]
    return p, u


def finite_rows(*arrays):
    """Number of leading samples in which every array is finite."""
    ok = np.ones(len(arrays[0]), dtype=bool)
    for a in arrays:
        ok &= np.isfinite(a.reshape(len(a), -1)).all(axis=1)
    return len(ok) if ok.all() else int(np.argmin(ok))


def growth(t, values):
    """Least-squares growth rate of log(values) over the last quarter."""
    keep = finite_rows(values)
    t, values = t[:keep], values[:keep]
    if keep < 4:
        return None, None, "zero"
    last = t >= t[0] + 0.75 * (t[-1] - t[0])
    v = values[last]
    if np.any(v <= 0) or last.sum() < 3:
        return None, None, "zero"
    rate = float(np.polyfit(t[last], np.log(v), 1)[0])
    factor = math.exp(rate * (t[last][-1] - t[last][0]))
    if factor > BOUND:
        kind = "growing"
    elif factor < 1.0 / BOUND:
        kind = "decaying"
    else:
        kind = "bounded"
    return rate, factor, kind


def main():
    s = json.loads((CASE / "caseSetup.json").read_text())
    dim, n, h, dt = s["dim"], s["n"], s["h"], s["deltaT"]
    half, a = s["domainHalfWidth"], s["halfWidth"]
    centres = -half + (np.arange(n) + 0.5) * h
    inside = np.abs(centres) < a
    m = int(inside.sum())

    # Cell volume and face area; 2D meshes are one cell of depth h thick.
    volume, area = h**3, h**2

    # Media from the initial alpha field. blockMesh orders the cells x
    # fastest, then y, then z, like the probe locations, so the arrays below
    # are indexed [z, y, x].
    alpha = np.zeros((n,) * dim)
    if s["dropRadius"] > 0:
        text = (CASE / "0" / "alpha.water").read_text()
        body = text.split("internalField", 1)[1]
        if "nonuniform" in body.split(";", 1)[0]:
            values = body.split("(", 1)[1].split(")", 1)[0].split()
            alpha = np.array(values, dtype=float).reshape((n,) * dim)
    rho_all = alpha * s["rhoL"] + (1 - alpha) * s["rhoG"]
    beta_all = (alpha / (s["rhoL"] * s["cL"] ** 2)
                + (1 - alpha) / (s["rhoG"] * s["cG"] ** 2))
    sel = np.ix_(*[np.where(inside)[0]] * dim)
    rho, beta = rho_all[sel].ravel(), beta_all[sel].ravel()

    t, p = read_probes("cellProbes", "p")
    _, q = read_probes("cellProbes", "pDot")
    _, u = read_probes("cellProbes", "U", components=dim)
    # A run that blew up is evaluated up to its last finite sample.
    keep = min(len(t), len(q), len(u))
    t, p, q, u = t[:keep], p[:keep], q[:keep], u[:keep]
    finite_samples = finite_rows(p, q, u)
    blowup = None if finite_samples == keep else float(t[finite_samples])
    t, p, q, u = (a[:finite_samples] for a in (t, p, q, u))
    steps = np.rint(t / dt).astype(int)
    if not np.array_equal(steps, s["probeStride"] * np.arange(1, len(t) + 1)):
        raise RuntimeError("probe times are not uniformly strided time steps")

    # Initial state at the cell centres.
    grids = np.meshgrid(*[centres[inside]] * dim, indexing="ij")
    # Transposing the (x, y[, z]) grid gives the probe order, x fastest.
    points = np.stack([g.T.ravel() for g in grids], axis=1)
    sq = np.sum((points - np.array(s["source"])[None, :dim]) ** 2, axis=1)
    p_init = s["p0"] * np.exp(-sq / (2 * s["pulseWidth"] ** 2))
    t = np.concatenate(([0.0], t))
    p = np.vstack((p_init, p))
    q = np.vstack((np.zeros(len(p_init)), q))
    u = np.concatenate((np.zeros((1, *u.shape[1:])), u))

    # Physical energy and discrete wave energy; the face term of E_h uses the
    # faces between two physical cells with the harmonic face density.
    e_ac = 0.5 * volume * (np.sum(beta * p**2, axis=1)
                           + np.sum(rho * np.sum(u**2, axis=2), axis=1))
    shape = (len(t),) + (m,) * dim
    p_grid, rho_grid = p.reshape(shape), rho.reshape((m,) * dim)
    stiffness = np.zeros(len(t))
    for axis in range(dim):
        lo = [slice(None)] * dim
        hi = [slice(None)] * dim
        lo[axis], hi[axis] = slice(None, -1), slice(1, None)
        inv_rho_f = 0.5 * (1 / rho_grid[tuple(lo)] + 1 / rho_grid[tuple(hi)])
        jump = p_grid[(slice(None), *hi)] - p_grid[(slice(None), *lo)]
        stiffness += np.sum((area / h) * inv_rho_f[None] * jump**2,
                            axis=tuple(range(1, dim + 1)))
    e_h = 0.5 * volume * np.sum(beta * q**2, axis=1) + 0.5 * stiffness

    summary = {"setup": s, "blowupTime": blowup,
               "lastSampleTime": float(t[-1])}
    obs = t <= T_OBS * (1 + 1e-9)
    if s["dropRadius"] == 0:
        p_ex, u_ex = free_space(s, points, t[obs])
        wp = (volume * beta)[None, :]
        e2 = {"pressure": math.sqrt(np.sum(wp * (p[obs] - p_ex) ** 2)
                                    / np.sum(wp * p_ex**2))}
        if u_ex is not None:
            wu = (volume * rho)[None, :, None]
            den = np.sum(wu * u_ex**2)
            e2["velocity"] = math.sqrt(np.sum(wu * (u[obs] - u_ex) ** 2) / den)
        e2["tObs"] = float(t[obs][-1])
        e2["samples"] = int(obs.sum())
        summary["E2exact"] = e2

    late = (t >= T_RESIDUAL) & (t <= T_OBS * (1 + 1e-9))
    summary["residual"] = {
        "from": T_RESIDUAL, "to": T_OBS,
        "maxAbsP": float(np.abs(p[late]).max() / s["p0"]) if late.any() else None,
        "maxAbsPAfter": float(np.abs(p[t >= T_RESIDUAL]).max() / s["p0"])
        if (t >= T_RESIDUAL).any() else None,
    }

    energy0 = {"Eac": float(e_ac[0]), "Eh": float(e_h[0])}
    after = t >= T_OBS
    summary["energy"] = {
        "initial": energy0,
        "atTobs": {"Eac": float(np.interp(T_OBS, t, e_ac)) / energy0["Eac"],
                   "Eh": float(np.interp(T_OBS, t, e_h)) / energy0["Eh"]}
        if t[-1] >= T_OBS else None,
        "maxAfterTobs": {"Eac": float(e_ac[after].max() / energy0["Eac"]),
                         "Eh": float(e_h[after].max() / energy0["Eh"])}
        if after.any() else None,
        "final": {"Eac": float(e_ac[-1] / energy0["Eac"]),
                  "Eh": float(e_h[-1] / energy0["Eh"])},
    }

    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "energy.csv").open("w") as stream:
        stream.write("t,Eac,Eh\n")
        for row in zip(t, e_ac, e_h):
            stream.write(",".join(f"{v:.10e}" for v in row) + "\n")

    classes = {"Eac": growth(t, e_ac), "Eh": growth(t, e_h)}
    domain = CASE / "postProcessing" / "domainProbes" / "0"
    if domain.is_dir():
        td, pd = read_probes("domainProbes", "p")
        _, vd = read_probes("domainProbes", "v")
        _, wd = read_probes("domainProbes", "u", components=dim)
        maxima = {"p": np.abs(pd).max(axis=1), "psi": np.abs(vd).max(axis=1),
                  "w": np.hypot.reduce(wd, axis=2).max(axis=1)}
        with (OUT / "maxima.csv").open("w") as stream:
            stream.write("t,maxAbsP,maxAbsPsi,maxAbsW\n")
            for i, ti in enumerate(td):
                stream.write(f"{ti:.10e}," + ",".join(
                    f"{maxima[k][i]:.6e}" for k in ("p", "psi", "w")) + "\n")
        summary["maxima"] = {
            k: {"overall": float(np.nanmax(v)), "final": float(v[-1])}
            for k, v in maxima.items()}
        for k, v in maxima.items():
            classes[f"max_{k}"] = growth(td, v)
    summary["lastQuarter"] = {
        k: {"rate": r, "factor": f, "class": c}
        for k, (r, f, c) in classes.items()}
    summary["finite"] = blowup is None

    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if blowup is not None:
        print(f"non-finite fields from t = {blowup:.4e} s")
    e2 = summary.get("E2exact", {})
    print(f"E2 exact: p {e2.get('pressure', float('nan')):.3e}, "
          f"U {e2.get('velocity', float('nan')):.3e}; "
          f"E_ac(end)/E_ac(0) {summary['energy']['final']['Eac']:.3e}, "
          f"E_h(end)/E_h(0) {summary['energy']['final']['Eh']:.3e}")
    r = summary["residual"]
    if r["maxAbsP"] is not None:
        print(f"residual max|p|/P0 for {T_RESIDUAL:g} <= t <= {T_OBS:g}: "
              f"{r['maxAbsP']:.3e}")
    for k, v in summary["lastQuarter"].items():
        if v["rate"] is not None:
            print(f"  {k:8s} {v['class']:9s} rate {v['rate']:.3e} 1/s, "
                  f"factor {v['factor']:.4g}")


if __name__ == "__main__":
    main()
