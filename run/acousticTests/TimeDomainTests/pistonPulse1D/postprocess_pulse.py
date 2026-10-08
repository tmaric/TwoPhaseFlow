#!/usr/bin/env python3
"""Piston-pulse verification for acousticWaveFoam.

Compares the probed pressure and velocity with the exact pulse solution
(pulseReference.py) over the physical domain x < pmlStart. E2 uses the
energy weights V*beta for pressure and V*rho for velocity, over the whole
physical domain and per medium. The PML is measured separately: the
largest pressure in the physical part of medium 2 at times when the exact
pressure vanishes there, relative to the largest transmitted pressure.
Writes postProcessing/pulseVerification/{summary.json,traces.csv}.
"""
import json
import math
from pathlib import Path

import numpy as np

from pulseReference import LayeredPulseSolution, PistonPulse

CASE = Path(__file__).resolve().parent
OUT = CASE / "postProcessing" / "pulseVerification"


def read_probes(field, nx):
    path = CASE / "postProcessing" / "cellProbes" / "0" / field
    times, rows = [], []
    with path.open() as stream:
        for line in stream:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.replace("(", " ").replace(")", " ").split()
            times.append(float(parts[0]))
            rows.append([float(v) for v in parts[1:]])
    data = np.array(rows)
    return np.array(times), (data if field == "p" else data[:, 0::3])


def e2(field, reference, weight):
    den = np.sum(weight * reference**2)
    return math.sqrt(np.sum(weight * (field - reference) ** 2) / den) if den else None


def main():
    s = json.loads((CASE / "caseSetup.json").read_text())
    nx, h, dt = s["nx"], s["h"], s["deltaT"]
    x = (np.arange(nx) + 0.5) * h

    t, p = read_probes("p", nx)
    _, u = read_probes("U", nx)
    if t[0] > 0.5 * dt:
        t = np.concatenate(([0.0], t))
        p = np.vstack((np.zeros(nx), p))
        u = np.vstack((np.zeros(nx), u))
    steps = np.rint(t / dt).astype(int)
    if not np.array_equal(steps, s["probeStride"] * np.arange(len(steps))):
        raise RuntimeError("probe times are not uniformly strided time steps")

    pulse = PistonPulse(s["f0"], s["nCycles"], s["u0"])
    exact = LayeredPulseSolution(pulse, s["rho1"], s["c1"], s["rho2"], s["c2"],
                                 s["xInterface"])
    p_ex, u_ex = exact.fields(x, t)

    medium1 = x < s["xInterface"]
    physical = x < s["pmlStart"]
    medium2 = physical & ~medium1
    rho = np.where(medium1, s["rho1"], s["rho2"])
    beta = 1.0 / np.where(medium1, s["rho1"] * s["c1"] ** 2,
                          s["rho2"] * s["c2"] ** 2)
    volume = h**3

    def weights(mask, coefficient):
        return np.where(mask, volume * coefficient, 0.0)[None, :]

    # Spurious pressure in medium 2 while the exact pressure vanishes there,
    # counted only once a transmitted pulse has passed.
    transmitted = float(np.abs(p_ex[:, medium2]).max())
    arrived = np.cumsum(np.abs(p_ex[:, medium2]).max(axis=1) > 0) > 0
    quiet = arrived & np.all(np.abs(p_ex[:, medium2]) <= 1e-9 * transmitted,
                             axis=1)
    spurious = (float(np.abs(p[quiet][:, medium2]).max() / transmitted)
                if transmitted > 0 and quiet.any() else None)

    summary = {
        "setup": s,
        "E2": {
            "pressure": e2(p, p_ex, weights(physical, beta)),
            "velocity": e2(u, u_ex, weights(physical, rho)),
            "pressureMedium1": e2(p, p_ex, weights(medium1, beta)),
            "pressureMedium2": e2(p, p_ex, weights(medium2, beta)),
            "velocityMedium1": e2(u, u_ex, weights(medium1, rho)),
            "velocityMedium2": e2(u, u_ex, weights(medium2, rho)),
        },
        "pml": {
            "transmittedAmplitude": float(transmitted),
            "quietSamples": int(quiet.sum()),
            "spuriousRelativePressure": spurious,
        },
        "maxAbsPressure": float(np.abs(p).max()),
        "finite": bool(np.isfinite(p).all() and np.isfinite(u).all()),
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    probes = [np.argmin(abs(x - xp)) for xp in
              (0.5 * s["xInterface"], s["xInterface"] - 2 * h,
               s["xInterface"] + 2 * h, 0.5 * (s["xInterface"] + s["pmlStart"]))]
    header = ["t"] + [f"{name}_x{x[i]:.5f}" for i in probes
                      for name in ("p", "pExact")]
    np.savetxt(OUT / "traces.csv",
               np.column_stack([t] + [col for i in probes
                                      for col in (p[:, i], p_ex[:, i])]),
               delimiter=",", header=",".join(header), comments="")

    e = {k: ("--" if v is None else f"{v:.3e}") for k, v in summary["E2"].items()}
    print(f"{s['config']}: nx = {nx}, N_lambda_g = {s['gasCellsPerWavelength']:.1f}, "
          f"N_T = {s['stepsPerPeriod']}, method {s['method']}")
    print(f"E2 pressure / velocity            {e['pressure']} {e['velocity']}")
    print(f"E2 pressure medium 1 / medium 2   {e['pressureMedium1']} "
          f"{e['pressureMedium2']}")
    print(f"E2 velocity medium 1 / medium 2   {e['velocityMedium1']} "
          f"{e['velocityMedium2']}")
    print(f"spurious pressure in quiet medium 2 {spurious} "
          f"({summary['pml']['quietSamples']} samples)")
    print(f"finite {summary['finite']}, max |p| {summary['maxAbsPressure']:.4e} Pa")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
