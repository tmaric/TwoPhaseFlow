#!/usr/bin/env python3
"""Standing-mode verification on the (non-orthogonal) square cavity.

Compares the probed p and U with the exact mode
  p = PA phi(x, y) cos(omega t),  U = -(PA/(rho omega)) grad phi sin(omega t),
  phi = cos(m pi x/L) cos(n pi y/L),
at the cell centres written by writeCellCentres, with volume weights:
E2 over all samples and over the last period. The modal amplitude
a(t) = sum V p phi / sum V phi^2 gives the realized amplification per
sample (fit of a[m+1] = c1 a[m] + c2 a[m-1]), hence a growth or decay rate
per period and the frequency error. The physical energy
E_ac = 1/2 sum V (beta p^2 + rho |U|^2), averaged over the first and the
last period, measures its drift.
Writes postProcessing/cavity2dVerification/summary.json.
"""
import json
import math
from pathlib import Path

import numpy as np

CASE = Path(__file__).resolve().parent
OUT = CASE / "postProcessing" / "cavity2dVerification"


def read_field(name, vector=False):
    text = (CASE / "0" / name).read_text().split("internalField", 1)[1].split(";", 1)[0]
    body = text.split("(", 1)[1].rsplit(")", 1)[0]
    data = np.array(body.replace("(", " ").replace(")", " ").split(), dtype=float)
    return data.reshape(-1, 3) if vector else data


def read_probes(field):
    times, rows = [], []
    with (CASE / "postProcessing/cellProbes/0" / field).open() as stream:
        for line in stream:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.replace("(", " ").replace(")", " ").split()
            times.append(float(parts[0]))
            rows.append(np.array(parts[1:], dtype=float))
    return np.array(times), np.array(rows)


def main():
    s = json.loads((CASE / "caseSetup.json").read_text())
    centres = read_field("C", vector=True)
    volume = read_field("V")
    x, y = centres[:, 0], centres[:, 1]
    length, m, n = s["L"], s["modeM"], s["modeN"]
    kx, ky = m * math.pi / length, n * math.pi / length
    omega, rho, c, pa = s["omega"], s["rho"], s["c"], s["Pa"]
    beta = 1 / (rho * c * c)
    phi = np.cos(kx * x) * np.cos(ky * y)
    grad = np.stack((-kx * np.sin(kx * x) * np.cos(ky * y),
                     -ky * np.cos(kx * x) * np.sin(ky * y)), axis=1)

    t, p = read_probes("p")
    _, u = read_probes("U")
    u = u.reshape(len(t), -1, 3)[:, :, :2]
    t = np.concatenate(([0.0], t))
    p = np.vstack((pa * phi, p))
    u = np.concatenate((np.zeros((1, *u.shape[1:])), u))
    finite = int(np.argmin(np.isfinite(p).all(axis=1))) if not np.isfinite(p).all() else len(t)
    t, p, u = t[:finite], p[:finite], u[:finite]

    p_ex = pa * np.outer(np.cos(omega * t), phi)
    u_ex = -(pa / (rho * omega)) * np.sin(omega * t)[:, None, None] * grad[None]

    def e2(a, b, w):
        w = w if a.ndim == 2 else w[..., None]
        return math.sqrt(np.sum(w * (a - b) ** 2) / np.sum(w * b ** 2))

    per = s["stepsPerPeriod"] // s["probeStride"]
    last = slice(len(t) - per, len(t))
    amp = (p * volume * phi).sum(axis=1) / (volume * phi * phi).sum()
    lhs = np.column_stack((amp[1:-1], amp[:-2]))
    (c1, c2), *_ = np.linalg.lstsq(lhs, amp[2:], rcond=None)
    modulus = math.sqrt(max(-c2, 0.0))
    arg = math.acos(max(-1.0, min(1.0, c1 / (2 * modulus)))) if modulus else 0.0
    samples_per_period = per
    growth_per_period = modulus ** samples_per_period
    omega_fit = arg / (s["probeStride"] * s["deltaT"])
    energy = 0.5 * (volume * (beta * p**2 + rho * (u**2).sum(axis=2))).sum(axis=1)
    first = slice(1, per + 1)
    summary = {
        "setup": s,
        "lastTime": float(t[-1]),
        "finite": finite == len(p) and bool(np.isfinite(p).all()),
        "E2": {"pressure": e2(p, p_ex, volume[None]),
               "velocity": e2(u, u_ex, volume[None]),
               "pressureLastPeriod": e2(p[last], p_ex[last], volume[None])},
        "modal": {"growthPerPeriod": growth_per_period,
                  "relativeFrequencyError": omega_fit / omega - 1.0,
                  "amplitudeFirst": float(np.abs(amp[first]).max()),
                  "amplitudeLast": float(np.abs(amp[last]).max())},
        "energyDrift": float(energy[last].mean() / energy[first].mean() - 1.0),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"E2 p {summary['E2']['pressure']:.3e}, U {summary['E2']['velocity']:.3e}, "
          f"last period p {summary['E2']['pressureLastPeriod']:.3e}")
    print(f"modal growth per period {growth_per_period - 1:+.3e}, "
          f"frequency error {summary['modal']['relativeFrequencyError']:+.3e}, "
          f"energy drift {summary['energyDrift']:+.3e}")


if __name__ == "__main__":
    main()
