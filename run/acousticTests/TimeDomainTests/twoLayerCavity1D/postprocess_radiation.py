#!/usr/bin/env python3
"""Radiation-average verification for the two-layer cavity mode.

Compares the probed running averages pr and momFlux_xx with three
references, averaged like the solver over the last N_s steps,
N_s = round(RADIATION_PERIODS/(f dt)), or the available steps 1..n:

  discrete   the fully discrete theta-rule solution from the initial
             pressure, expanded in all discrete modes of the assembled
             operator, with U from the theta velocity update of the
             face-flux reconstruction;
  exact      the exact mode sampled at every step;
  mean       the exact time means, pr = beta P^2/4 - rho Uhat^2/4 and
             R_xx = rho Uhat^2/2.

Errors are maxima over the sampled times with a full window and over four
cell groups: gas and liquid interiors and the gas and liquid cells next to
the interface. Each is relative to the largest exact mean of that medium.
Serial runs with the interface on a face only.
Writes postProcessing/radiationVerification/summary.json.
"""
import json
import math
from pathlib import Path

import numpy as np

from postprocess_two_layer import internal_field, latest_time, read_probes
from twoLayerModes import TwoLayerMode

CASE = Path(__file__).resolve().parent
OUT = CASE / "postProcessing" / "radiationVerification"


def window_average(values, n_s, rows):
    """Running mean of values[1..n] over the last n_s steps at steps `rows`."""
    csum = np.vstack((np.zeros_like(values[:1]),
                      np.cumsum(values[1:], axis=0)))
    lo = np.maximum(rows - n_s, 0)
    count = np.maximum(rows - lo, 1)[:, None]
    return (csum[rows] - csum[lo]) / count


def main():
    s = json.loads((CASE / "caseSetup.json").read_text())
    if s["nProcs"] != 1 or s["interfaceOffset"] != 0:
        raise SystemExit("radiation check: serial, interface on a face")
    nx, h, dt, theta = s["nx"], s["h"], s["deltaT"], s["theta"]
    n_gas, n_steps = s["nGas"], s["nSteps"]
    n_s = max(1, round(s["radiationPeriods"] / (s["f"] * dt)))
    x = (np.arange(nx) + 0.5) * h

    alpha = internal_field(CASE / "0" / "alpha.water", nx)
    rho = alpha * s["rhoL"] + (1 - alpha) * s["rhoG"]
    beta = (alpha / (s["rhoL"] * s["cL"] ** 2)
            + (1 - alpha) / (s["rhoG"] * s["cG"] ** 2))
    a = internal_field(latest_time(CASE) / "invRhof", nx - 1)

    # beta q' = -(1/h^2) K p with the solver's face coefficients.
    k = np.zeros((nx, nx))
    for f in range(nx - 1):
        k[f, f] += a[f]
        k[f, f + 1] -= a[f]
        k[f + 1, f + 1] += a[f]
        k[f + 1, f] -= a[f]
    scale = 1.0 / np.sqrt(beta)
    lam, w = np.linalg.eigh(scale[:, None] * k * scale[None, :] / h**2)
    xi = np.sqrt(np.clip(lam, 0.0, None)) * dt
    g = (1 + (1 - theta) * 1j * xi) / (1 - theta * 1j * xi)
    p0 = internal_field(CASE / "0" / "p", nx)
    coeff = w.T @ (p0 / scale)
    steps = np.arange(n_steps + 1)
    amp = np.real(np.power.outer(g, steps)).T * coeff[None, :]
    p_d = amp @ (w.T * scale[None, :])

    # Face-flux reconstruction in pure cells:
    # [a_e (p_E - p_P) + a_w (p_P - p_W)]/(2 h), zero flux at the walls.
    recon = np.zeros((nx, nx))
    for f in range(nx - 1):
        for cell, sign in ((f, 1.0), (f + 1, 1.0)):
            recon[cell, f + 1] += sign * a[f] / (2 * h)
            recon[cell, f] -= sign * a[f] / (2 * h)
    accel = p_d @ recon.T
    increments = theta * accel[1:] + (1 - theta) * accel[:-1]
    u_d = -dt * np.vstack((np.zeros(nx), np.cumsum(increments, axis=0)))

    mode = TwoLayerMode(s["rhoG"], s["cG"], s["rhoL"], s["cL"],
                        s["xInterface"], s["length"], s["omega"], s["Pa"])
    t_all = steps * dt
    shape_p, shape_u = mode.pressure(x), mode.velocity(x)
    p_e = np.outer(np.cos(s["omega"] * t_all), shape_p)
    u_e = np.outer(np.sin(s["omega"] * t_all), shape_u)

    t, pr = read_probes("pr")
    _, mom = read_probes("momFlux")
    _, p_probe = read_probes("p")
    _, u_probe = read_probes("U")
    rows = np.rint(t / dt).astype(int)
    r_xx = mom[:, 0::9]
    u_probe = u_probe[:, 0::3]
    full = rows >= n_s

    def averages(p, u):
        return (window_average(0.5 * beta * p**2 - 0.5 * rho * u**2, n_s, rows),
                window_average(rho * u**2, n_s, rows))

    mean_pr = 0.25 * beta * shape_p**2 - 0.25 * rho * shape_u**2
    mean_r = 0.5 * rho * shape_u**2
    refs = {"discrete": averages(p_d, u_d), "exact": averages(p_e, u_e),
            "mean": (np.broadcast_to(mean_pr, pr.shape),
                     np.broadcast_to(mean_r, pr.shape))}

    gas = np.arange(nx) < n_gas
    groups = {
        "gasInterior": np.arange(nx) < n_gas - 1,
        "gasInterfaceCell": np.arange(nx) == n_gas - 1,
        "liquidInterfaceCell": np.arange(nx) == n_gas,
        "liquidInterior": np.arange(nx) > n_gas,
    }
    scales = {}
    for medium, mask in (("gas", gas), ("liquid", ~gas)):
        scales[medium] = (np.abs(mean_pr[mask]).max(), mean_r[mask].max())

    errors = {}
    for name, (ref_pr, ref_r) in refs.items():
        errors[name] = {}
        for group, mask in groups.items():
            sp, sr = scales["gas" if group.startswith("gas") else "liquid"]
            errors[name][group] = {
                "pr": float(np.abs(pr - ref_pr)[full][:, mask].max() / sp),
                "Rxx": float(np.abs(r_xx - ref_r)[full][:, mask].max() / sr),
            }
    fields = {
        "p": float(np.abs(p_probe - p_d[rows]).max() / np.abs(p0).max()),
        "U": float(np.abs(u_probe - u_d[rows]).max()
                   / np.abs(u_d).max()),
    }

    summary = {"setup": s, "windowSamples": n_s, "errors": errors,
               "fieldsVsDiscrete": fields,
               "meanScales": {m: list(v) for m, v in scales.items()}}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"mode {s['mode']}, nx {nx}, N_T {s['stepsPerPeriod']}, "
          f"N_s {n_s}; fields vs discrete: p {fields['p']:.3e}, "
          f"U {fields['U']:.3e}")
    for name, by_group in errors.items():
        print(f"  {name}")
        for group, e in by_group.items():
            print(f"    {group:20s} pr {e['pr']:.3e}  Rxx {e['Rxx']:.3e}")


if __name__ == "__main__":
    main()
