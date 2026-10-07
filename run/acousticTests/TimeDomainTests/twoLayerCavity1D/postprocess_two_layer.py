#!/usr/bin/env python3
"""Two-layer cavity verification for acousticWaveFoam.

The discrete reference is built from the face coefficients the solver wrote
(invRhof) and the cell coefficients of alpha.water, so it represents the
operator actually solved, whatever the interface closure. Each processor
side of a decomposed face contributes to its own cell row; unequal sides
make the operator nonsymmetric and are reported.

The probed pressure is projected onto the discrete eigenmode nearest to the
exact frequency with the left eigenvector, so other modes in the initial
field do not affect the projected amplitude. Writes
postProcessing/twoLayerVerification/{summary.json,history.csv}.
"""
import json
import math
import re
from pathlib import Path

import numpy as np

from twoLayerModes import TwoLayerMode

CASE = Path(__file__).resolve().parent
OUT = CASE / "postProcessing" / "twoLayerVerification"


def read_probes(field):
    path = CASE / "postProcessing" / "cellProbes" / "0" / field
    times, rows = [], []
    with path.open() as stream:
        for line in stream:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.replace("(", " ").replace(")", " ").split()
            times.append(float(parts[0]))
            rows.append([float(v) for v in parts[1:]])
    return np.array(times), np.array(rows)


def internal_field(path, size):
    text = path.read_text()
    match = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*"
                      r"\((.*?)\)\s*;", text, re.S)
    if match:
        return np.array(match.group(2).split(), dtype=float)
    match = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", text)
    return np.full(size, float(match.group(1)))


def patch_value(path, patch):
    text = path.read_text()
    block = re.search(re.escape(patch) + r"\s*\{(.*?)\}", text, re.S).group(1)
    match = re.search(r"value\s+nonuniform\s+List<scalar>\s*\d+\s*\((.*?)\)",
                      block, re.S)
    if match:
        return np.array(match.group(1).split(), dtype=float)
    match = re.search(r"value\s+uniform\s+([-+0-9.eE]+)", block)
    return np.array([float(match.group(1))])


def latest_time(root):
    times = [d for d in root.iterdir() if d.is_dir()
             and re.fullmatch(r"[0-9.eE+-]+", d.name) and float(d.name) > 0]
    return max(times, key=lambda d: float(d.name))


def fit_amplification(a, stride):
    lhs = np.column_stack((a[1:-1], a[:-2]))
    (c1, c2), *_ = np.linalg.lstsq(lhs, a[2:], rcond=None)
    modulus = math.sqrt(max(-c2, 0.0))
    argument = math.acos(max(-1.0, min(1.0, c1 / (2.0 * modulus))))
    modulus, argument = modulus ** (1.0 / stride), argument / stride
    return modulus * complex(math.cos(argument), math.sin(argument))


def discrete_reference(method, theta, xi, steps):
    if method == "fullyImplicit":
        a = np.empty(steps[-1] + 1)
        a[0] = 1.0
        previous = 1.0
        for m in range(steps[-1]):
            a[m + 1] = (2.0 * a[m] - previous) / (1.0 + xi * xi)
            previous = a[m]
        g = (1.0 + 1j * xi) / (1.0 + xi * xi)
        return a[steps], g
    g = (1.0 + (1.0 - theta) * 1j * xi) / (1.0 - theta * 1j * xi)
    return np.real(g ** steps), g


def main():
    s = json.loads((CASE / "caseSetup.json").read_text())
    nx, h, area = s["nx"], s["h"], s["crossSection"]
    volume = area * h
    method, theta, dt = s["method"], s["theta"], s["deltaT"]
    x = (np.arange(nx) + 0.5) * h

    alpha = internal_field(CASE / "0" / "alpha.water", nx)
    rho = alpha * s["rhoL"] + (1 - alpha) * s["rhoG"]
    beta = (alpha / (s["rhoL"] * s["cL"] ** 2)
            + (1 - alpha) / (s["rhoG"] * s["cG"] ** 2))

    # Face coefficients seen by the owner row and by the neighbour row.
    t_end = latest_time(CASE)
    a_face = internal_field(t_end / "invRhof", nx - 1)
    a_own, a_nei = a_face.copy(), a_face.copy()
    proc_mismatch = None
    if s["nProcs"] == 2:
        n_gas = s["nGas"]
        t0 = CASE / "processor0" / t_end.name
        t1 = CASE / "processor1" / t_end.name
        side0 = patch_value(t0 / "invRhof", "procBoundary0to1")[0]
        side1 = patch_value(t1 / "invRhof", "procBoundary1to0")[0]
        alpha0 = patch_value(t0 / "alphaf", "procBoundary0to1")[0]
        alpha1 = patch_value(t1 / "alphaf", "procBoundary1to0")[0]
        face = n_gas - 1
        a_own[face], a_nei[face] = side0, side1
        proc_mismatch = {
            "invRhofRank0": side0, "invRhofRank1": side1,
            "alphafRank0": alpha0, "alphafRank1": alpha1,
            "relativeMismatch": abs(side0 - side1) / max(side0, side1),
        }

    # beta p'' = -(1/h^2) K p, rows assembled side by side.
    k = np.zeros((nx, nx))
    for f in range(nx - 1):
        k[f, f] += a_own[f]
        k[f, f + 1] -= a_own[f]
        k[f + 1, f + 1] += a_nei[f]
        k[f + 1, f] -= a_nei[f]
    omega = s["omega"]
    if np.array_equal(a_own, a_nei):
        # Symmetric pencil: B^{-1/2} K B^{-1/2} is symmetric, so eigh keeps
        # the B-orthogonality of the modes to round-off.
        scale = 1.0 / np.sqrt(beta)
        eigenvalues, vectors = np.linalg.eigh(
            scale[:, None] * k * scale[None, :] / (h * h))
        index = np.argmin(abs(np.sqrt(np.abs(eigenvalues)) - omega))
        omega_h = math.sqrt(eigenvalues[index])
        shape = scale * vectors[:, index]
        projector = vectors[:, index] / scale
    else:
        eigenvalues, right = np.linalg.eig(k / (h * h * beta[:, None]))
        left = np.linalg.inv(right)
        index = np.argmin(abs(np.sqrt(np.abs(eigenvalues.real)) - omega))
        omega_h = math.sqrt(eigenvalues[index].real)
        projector = left[index].real
        shape = right[:, index].real

    t, p = read_probes("p")
    _, q = read_probes("pDot")
    _, u = read_probes("U")
    u = u.reshape(len(t), nx, 3)[:, :, 0]
    if t[0] > 0.5 * dt:
        t = np.concatenate(([0.0], t))
        p = np.vstack((internal_field(CASE / "0" / "p", nx), p))
        q = np.vstack((np.zeros(nx), q))
        u = np.vstack((np.zeros(nx), u))
    stride = s["probeStride"]
    steps = np.rint(t / dt).astype(int)
    if not np.array_equal(steps, stride * np.arange(len(steps))):
        raise RuntimeError("probe times are not uniformly strided time steps")
    if method == "fullyImplicit":
        q = (np.vstack((np.zeros(nx), np.diff(p, axis=0) / dt))
             if stride == 1 else None)

    amplitude = p @ projector
    a0 = amplitude[0]
    amplitude = amplitude / a0
    mass = volume * beta
    modal_initial = a0 * shape
    modal_share = math.sqrt(
        np.sum(mass * modal_initial**2) / np.sum(mass * p[0] ** 2))

    xi = omega_h * dt
    a_disc, g_disc = discrete_reference(method, theta, xi, steps)
    g_fit = fit_amplification(amplitude, stride)
    omega_fit = math.atan2(g_fit.imag, g_fit.real) / dt
    gamma_fit = -math.log(abs(g_fit)) / dt
    omega_h_fit = (2.0 / dt * math.tan(0.5 * omega_fit * dt)
                   if theta == 0.5 and method != "fullyImplicit" else None)

    mode = TwoLayerMode(s["rhoG"], s["cG"], s["rhoL"], s["cL"],
                        s["xInterface"], s["length"], omega, s["Pa"])
    p_exact = np.outer(np.cos(omega * t), mode.pressure(x))
    u_exact = np.outer(np.sin(omega * t), mode.velocity(x))

    def e2(field, reference, weight):
        return math.sqrt(np.sum(weight * (field - reference) ** 2)
                         / np.sum(weight * reference**2))

    # Discrete wave energy with the face-averaged coefficient.
    a_mean = 0.5 * (a_own + a_nei)
    potential = 0.5 * area / h * np.sum(a_mean * np.diff(p, axis=1) ** 2,
                                        axis=1)
    energy_h = (0.5 * volume * np.sum(beta * q**2, axis=1) + potential
                if q is not None else np.full(len(t), np.nan))
    energy_ac = 0.5 * volume * np.sum(beta * p**2 + rho * u**2, axis=1)

    symmetric = bool(np.allclose(a_own, a_nei, rtol=0, atol=0))
    summary = {
        "setup": s,
        "interfaceFace": {
            "index": int(s["nGas"] - 1),
            "invRhofOwnerRow": float(a_own[s["nGas"] - 1]),
            "invRhofNeighbourRow": float(a_nei[s["nGas"] - 1]),
            "linearClosure": 2.0 / (s["rhoG"] + s["rhoL"]),
        },
        "processorFace": proc_mismatch,
        "operatorSymmetric": symmetric,
        "omega": omega,
        "omega_h": omega_h,
        "xi": xi,
        "initialModalShare": modal_share,
        "implementation": {
            "maxModalErrorVsDiscrete": float(np.abs(amplitude - a_disc).max()),
        },
        "rates": {
            "omegaFit": omega_fit,
            "omegaDiscrete": math.atan2(g_disc.imag, g_disc.real) / dt,
            "gammaFit": gamma_fit,
            "gammaDiscrete": -math.log(abs(g_disc)) / dt,
            "omegaHFromTrapezoidalInverse": omega_h_fit,
            "relativeFrequencyErrorDiscreteVsExact": omega_h / omega - 1.0,
            "relativeFrequencyErrorFitVsExact": (
                omega_h_fit / omega - 1.0 if omega_h_fit else None),
        },
        "E2": {
            "pressureVsExact": e2(p, p_exact, volume * beta),
            "velocityVsExact": e2(u, u_exact, volume * rho),
        },
        "energy": {
            "discreteWaveEnergyMaxRelativeDrift": (
                float(np.abs(energy_h / energy_h[0] - 1.0).max())
                if q is not None else None),
            "discreteWaveEnergyFinalRatio": (
                float(energy_h[-1] / energy_h[0]) if q is not None else None),
        },
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    np.savetxt(
        OUT / "history.csv",
        np.column_stack((t, amplitude, a_disc, np.cos(omega * t), energy_h,
                         energy_ac)),
        delimiter=",",
        header="t,a,aDiscrete,aExact,discreteWaveEnergy,physicalEnergy",
        comments="",
    )

    print(f"mode {s['mode']} f = {s['f']:.6f} Hz, nx = {nx}, "
          f"method {method}, theta {theta}, nProcs {s['nProcs']}")
    print(f"interface face invRhof own/nei  "
          f"{summary['interfaceFace']['invRhofOwnerRow']:.6g} / "
          f"{summary['interfaceFace']['invRhofNeighbourRow']:.6g}  "
          f"(linear closure {summary['interfaceFace']['linearClosure']:.6g})")
    if proc_mismatch:
        print(f"processor face invRhof rank0/rank1  "
              f"{proc_mismatch['invRhofRank0']:.6g} / "
              f"{proc_mismatch['invRhofRank1']:.6g}")
    print(f"operator symmetric          {symmetric}")
    print(f"initial modal share         {modal_share:.12f}")
    print(f"maxModalErrorVsDiscrete     "
          f"{summary['implementation']['maxModalErrorVsDiscrete']:.3e}")
    print(f"omega_h/omega - 1 (discrete) {omega_h / omega - 1.0:.6e}")
    if omega_h_fit:
        print(f"omega_h,fit/omega_h - 1      {omega_h_fit / omega_h - 1.0:.3e}")
    print(f"gamma fit / discrete        {gamma_fit:.6e} "
          f"{summary['rates']['gammaDiscrete']:.6e}")
    print(f"E2 pressure / velocity      {summary['E2']['pressureVsExact']:.3e} "
          f"{summary['E2']['velocityVsExact']:.3e}")
    print(f"discrete energy max drift   "
          f"{summary['energy']['discreteWaveEnergyMaxRelativeDrift']}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
