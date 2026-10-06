#!/usr/bin/env python3
"""Closed-cavity standing-mode verification for acousticWaveFoam.

On the uniform cell-centred mesh with rigid walls, cos(k x_P) is an exact
eigenvector of the orthogonal pressure operator with frequency
omega_h = (2 c / h) sin(k h / 2). The probed pressure is projected onto this
mode and compared with three references:

  exact            cos(omega t)
  semi-discrete    cos(omega_h t)
  fully discrete   the recurrence of the selected time integrator

Writes postProcessing/cavityVerification/{summary.json,history.csv}.
"""
import json
import math
import re
from pathlib import Path

import numpy as np

CASE = Path(__file__).resolve().parent
OUT = CASE / "postProcessing" / "cavityVerification"


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


def read_initial_pressure(nx):
    text = (CASE / "0" / "p").read_text()
    match = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\((.*?)\)\s*;",
        text,
        re.S,
    )
    if not match or int(match.group(1)) != nx:
        raise RuntimeError("0/p does not contain the nonuniform initial mode")
    return np.array(match.group(2).split(), dtype=float)


def fit_amplification(a, stride):
    """Fit a[m+1] = c1 a[m] + c2 a[m-1] on samples `stride` steps apart.

    The fitted factor is G**stride; the single-step G is unique while
    stride*arg(G) < pi, i.e. for more than two samples per period.
    """
    lhs = np.column_stack((a[1:-1], a[:-2]))
    (c1, c2), *_ = np.linalg.lstsq(lhs, a[2:], rcond=None)
    modulus = math.sqrt(max(-c2, 0.0))
    argument = math.acos(max(-1.0, min(1.0, c1 / (2.0 * modulus))))
    modulus, argument = modulus ** (1.0 / stride), argument / stride
    return modulus * complex(math.cos(argument), math.sin(argument))


def discrete_reference(method, theta, xi, steps):
    if method == "fullyImplicit":
        # Eq. (backwardPressure) with a^{-1} = a^0, i.e. zero initial rate.
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
    rho, c, nx, h, k = s["rho"], s["c"], s["nx"], s["h"], s["k"]
    pa, dt, method, theta = s["Pa"], s["deltaT"], s["method"], s["theta"]
    area = s["crossSection"]
    volume = area * h
    beta = 1.0 / (rho * c * c)

    x = (np.arange(nx) + 0.5) * h
    mode = np.cos(k * x)
    omega = c * k
    omega_h = 2.0 * c / h * math.sin(0.5 * k * h)
    xi = omega_h * dt

    t, p = read_probes("p")
    _, q = read_probes("pDot")
    _, u = read_probes("U")
    u = u.reshape(len(t), nx, 3)[:, :, 0]
    if p.shape[1] != nx:
        raise RuntimeError(f"expected {nx} probes, found {p.shape[1]}")

    if t[0] > 0.5 * dt:
        t = np.concatenate(([0.0], t))
        p = np.vstack((read_initial_pressure(nx), p))
        q = np.vstack((np.zeros(nx), q))
        u = np.vstack((np.zeros(nx), u))

    stride = s["probeStride"]
    steps = np.rint(t / dt).astype(int)
    if not np.array_equal(steps, stride * np.arange(len(steps))):
        raise RuntimeError("probe times are not uniformly strided time steps")

    if method == "fullyImplicit":
        # The backward branch does not update pDot; its own rate b^n is
        # the backward difference, available only when every step is probed.
        q = (
            np.vstack((np.zeros(nx), np.diff(p, axis=0) / dt))
            if stride == 1
            else None
        )

    amplitude = p @ mode / (mode @ mode) / pa
    modal_part = pa * np.outer(amplitude, mode)
    non_modal = np.linalg.norm(p - modal_part, axis=1).max() / (
        pa * np.linalg.norm(mode)
    )

    a_exact = np.cos(omega * t)
    a_semi = np.cos(omega_h * t)
    a_disc, g_disc = discrete_reference(method, theta, xi, steps)

    def e2(reference):
        field = pa * np.outer(reference, mode)
        return math.sqrt(
            np.sum(volume * (p - field) ** 2) / np.sum(volume * field**2)
        )

    u_exact = pa / (rho * c) * np.outer(np.sin(omega * t), np.sin(k * x))
    e2_velocity = math.sqrt(
        np.sum(volume * (u - u_exact) ** 2) / np.sum(volume * u_exact**2)
    )

    g_fit = fit_amplification(amplitude, stride)
    omega_fit = math.atan2(g_fit.imag, g_fit.real) / dt
    gamma_fit = -math.log(abs(g_fit)) / dt
    omega_disc = math.atan2(g_disc.imag, g_disc.real) / dt
    gamma_disc = -math.log(abs(g_disc)) / dt
    omega_h_fit = (
        2.0 / dt * math.tan(0.5 * omega_fit * dt) if theta == 0.5 else None
    )

    # Eq. (discreteWaveEnergy) with the orthogonal face coefficients A/(rho h).
    energy_h = (
        0.5 * volume * beta * np.sum(q**2, axis=1)
        + 0.5 * area / (rho * h) * np.sum(np.diff(p, axis=1) ** 2, axis=1)
        if q is not None
        else np.full(len(t), np.nan)
    )
    # Eq. (physicalEnergy).
    energy_ac = 0.5 * volume * np.sum(beta * p**2 + rho * u**2, axis=1)
    last_period = t >= t[-1] - 1.0 / s["f"]

    n_steps = int(steps[-1])
    summary = {
        "setup": s,
        "omega": omega,
        "omega_h": omega_h,
        "xi": xi,
        "courant": c * dt / h,
        "nonModalContentMax": non_modal,
        "implementation": {
            "maxModalErrorVsDiscrete": float(np.abs(amplitude - a_disc).max()),
            "maxFieldErrorVsDiscrete": float(
                np.abs(p - pa * np.outer(a_disc, mode)).max() / pa
            ),
        },
        "rates": {
            "omegaFit": omega_fit,
            "omegaDiscrete": omega_disc,
            "gammaFit": gamma_fit,
            "gammaDiscrete": gamma_disc,
            "omegaHFromTrapezoidalInverse": omega_h_fit,
            "relativeFrequencyErrorVsExact": omega_fit / omega - 1.0,
            "relativeFrequencyErrorVsSemiDiscrete": omega_fit / omega_h - 1.0,
            "relativeSpatialFrequencyError": (
                omega_h_fit / omega - 1.0 if omega_h_fit else None
            ),
            "relativeSpatialFrequencyErrorPredicted": omega_h / omega - 1.0,
        },
        "endOfRun": {
            "amplitudeFactorFit": abs(g_fit) ** n_steps,
            "amplitudeFactorDiscrete": abs(g_disc) ** n_steps,
            "phaseErrorVsExact_rad": (omega_fit - omega) * t[-1],
            "phaseErrorVsSemiDiscrete_rad": (omega_fit - omega_h) * t[-1],
        },
        "E2": {
            "vsExact": e2(a_exact),
            "vsSemiDiscrete": e2(a_semi),
            "vsDiscrete": e2(a_disc),
            "velocityVsExact": e2_velocity,
        },
        "energy": {
            "discreteWaveEnergyMaxRelativeDrift": (
                float(np.abs(energy_h / energy_h[0] - 1.0).max())
                if q is not None
                else None
            ),
            "discreteWaveEnergyFinalRatio": (
                float(energy_h[-1] / energy_h[0]) if q is not None else None
            ),
            "physicalEnergyFinalMeanRatio": float(
                energy_ac[last_period].mean() / energy_ac[0]
            ),
            "physicalEnergyLastPeriodRelativeSwing": float(
                np.ptp(energy_ac[last_period]) / energy_ac[last_period].mean()
            ),
        },
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    np.savetxt(
        OUT / "history.csv",
        np.column_stack(
            (t, amplitude, a_disc, a_semi, a_exact, energy_h, energy_ac)
        ),
        delimiter=",",
        header="t,a,aDiscrete,aSemiDiscrete,aExact,discreteWaveEnergy,"
        "physicalEnergy",
        comments="",
    )

    print(f"method {method}, theta {theta}, Courant {summary['courant']:.4g}")
    print(f"xi = omega_h dt = {xi:.6g}")
    print(f"non-modal content          {non_modal:.3e}")
    for key, value in summary["implementation"].items():
        print(f"{key:26s} {value:.3e}")
    print(f"E2 vs exact / semi / disc  {summary['E2']['vsExact']:.3e} "
          f"{summary['E2']['vsSemiDiscrete']:.3e} "
          f"{summary['E2']['vsDiscrete']:.3e}")
    print(f"E2 velocity vs exact       {e2_velocity:.3e}")
    print(f"gamma fit / discrete       {gamma_fit:.6e} {gamma_disc:.6e}")
    print(f"omega fit / discrete       {omega_fit:.10e} {omega_disc:.10e}")
    print(f"discrete energy max drift  "
          f"{summary['energy']['discreteWaveEnergyMaxRelativeDrift']}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
