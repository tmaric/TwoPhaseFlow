#!/usr/bin/env python3
"""Radiation-average verification for the closed-cavity mode.

Compares the probed running averages pr and momFlux_xx with three
references, each averaged like the solver: over the last N_s steps,
N_s = round(RADIATION_PERIODS/(f dt)), and over the available steps
1..n during startup (the initial state is not part of the window).

  discrete   the fully discrete mode of the theta rule: p = PA Re(g^n)
             cos(k x_P), with U from the theta velocity update applied to
             the discrete pressure gradient of the mode;
  exact      the exact mode sampled at the same times;
  mean       the exact time means, pr = beta PA^2/4 cos(2 k x) and
             R_xx = beta PA^2/2 sin^2(k x).

Errors are maxima over all cells and over the samples with a full window,
relative to beta PA^2/4 (pr) and beta PA^2/2 (R_xx). Also reported is the
bound |sin(N_s phi)|/(N_s sin phi) of the window bias, phi = arg(g).
Writes postProcessing/radiationVerification/summary.json.
"""
import json
import math
from pathlib import Path

import numpy as np

CASE = Path(__file__).resolve().parent
OUT = CASE / "postProcessing" / "radiationVerification"


def read_probes(field):
    path = CASE / "postProcessing" / "cellProbes" / "0" / field
    times, rows = [], []
    with path.open() as stream:
        for line in stream:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.replace("(", " ").replace(")", " ").split()
            times.append(float(parts[0]))
            rows.append(np.array(parts[1:], dtype=float))
    return np.array(times), np.array(rows)


def window_average(values, n_s):
    """Running mean of values[1..n] over the last n_s rows (row 0 = t = 0)."""
    csum = np.cumsum(values[1:], axis=0)
    csum = np.vstack((np.zeros_like(values[:1]), csum))
    n = np.arange(len(values))
    lo = np.maximum(n - n_s, 0)
    count = np.maximum(n - lo, 1)[:, None]
    return (csum - csum[lo]) / count


def main():
    s = json.loads((CASE / "caseSetup.json").read_text())
    if s["probeStride"] != 1:
        raise SystemExit("radiation check needs every step probed")
    if s["method"] == "fullyImplicit":
        raise SystemExit("radiation check covers the theta rule only")
    nx, h, rho, c = s["nx"], s["h"], s["rho"], s["c"]
    k, pa, dt, f, theta = s["k"], s["Pa"], s["deltaT"], s["f"], s["theta"]
    beta = 1.0 / (rho * c * c)
    n_s = max(1, round(s["radiationPeriods"] / (f * dt)))

    t, pr = read_probes("pr")
    _, mom = read_probes("momFlux")
    r_xx = mom[:, 0::9]
    n_steps = len(t)
    if not np.allclose(t, dt * np.arange(1, n_steps + 1), rtol=0, atol=1e-6 * dt):
        raise RuntimeError("probe times are not every time step")
    t = np.concatenate(([0.0], t))
    pr = np.vstack((np.zeros(nx), pr))
    r_xx = np.vstack((np.zeros(nx), r_xx))

    x = (np.arange(nx) + 0.5) * h
    mode = np.cos(k * x)

    # Fully discrete theta-rule mode.
    omega_h = 2 * c / h * math.sin(k * h / 2)
    xi = omega_h * dt
    g = (1 + (1 - theta) * 1j * xi) / (1 - theta * 1j * xi)
    amp = np.real(g ** np.arange(n_steps + 1))
    # Face-flux reconstruction on the uniform 1D mesh: (p_E - p_W)/(2 h rho),
    # with the zero-flux wall face standing in for the missing neighbour.
    east = np.append(mode[1:], mode[-1])
    west = np.insert(mode[:-1], 0, mode[0])
    grad_mode = (east - west) / (2 * h * rho)
    increments = theta * amp[1:] + (1 - theta) * amp[:-1]
    u_amp = -dt * np.concatenate(([0.0], np.cumsum(increments)))
    p_d = pa * np.outer(amp, mode)
    u_d = pa * np.outer(u_amp, grad_mode)

    # Exact mode.
    omega = c * k
    p_e = pa * np.outer(np.cos(omega * t), mode)
    u_e = pa / (rho * c) * np.outer(np.sin(omega * t), np.sin(k * x))

    def averages(p, u):
        return (window_average(0.5 * beta * p**2 - 0.5 * rho * u**2, n_s),
                window_average(rho * u**2, n_s))

    refs = {"discrete": averages(p_d, u_d), "exact": averages(p_e, u_e)}
    full = np.arange(n_steps + 1) >= n_s
    mean = (np.broadcast_to(0.25 * beta * pa**2 * np.cos(2 * k * x), pr.shape),
            np.broadcast_to(0.5 * beta * pa**2 * np.sin(k * x) ** 2, pr.shape))
    refs["mean"] = mean
    scale = (0.25 * beta * pa**2, 0.5 * beta * pa**2)

    errors = {}
    for name, (ref_pr, ref_r) in refs.items():
        errors[name] = {
            "pr": float(np.abs(pr - ref_pr)[full].max() / scale[0]),
            "Rxx": float(np.abs(r_xx - ref_r)[full].max() / scale[1]),
        }
    # The averaging alone: the solver's own p and U through the window.
    _, p_probe = read_probes("p")
    _, u_probe = read_probes("U")
    p_probe = np.vstack((pa * mode, p_probe))
    u_probe = np.vstack((np.zeros(nx), u_probe[:, 0::3]))
    own_pr, own_r = averages(p_probe, u_probe)
    errors["ownFields"] = {
        "pr": float(np.abs(pr - own_pr)[1:].max() / scale[0]),
        "Rxx": float(np.abs(r_xx - own_r)[1:].max() / scale[1]),
    }
    # Fields: discrete reference against the probed fields.
    fields = {
        "p": float(np.abs(p_probe - p_d).max() / pa),
        "U": float(np.abs(u_probe - u_d).max() / (pa / (rho * c))),
    }

    phi = abs(np.angle(g))
    bias_bound = abs(math.sin(n_s * phi)) / (n_s * math.sin(phi))
    exact_bound = abs(math.sin(n_s * omega * dt)) / (n_s * math.sin(omega * dt))
    summary = {
        "setup": s, "windowSamples": n_s,
        "windowDiscretePeriods": n_s * phi / (2 * math.pi),
        "biasBoundDiscrete": bias_bound, "biasBoundExact": exact_bound,
        "errors": errors, "fieldsVsDiscrete": fields,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"N_s = {n_s} ({summary['windowDiscretePeriods']:.6f} discrete "
          f"periods), bias bound {bias_bound:.3e} (exact-mode {exact_bound:.3e})")
    print(f"fields vs discrete: p {fields['p']:.3e}, U {fields['U']:.3e}")
    for name, e in errors.items():
        print(f"  {name:9s} pr {e['pr']:.3e}  Rxx {e['Rxx']:.3e}")


if __name__ == "__main__":
    main()
