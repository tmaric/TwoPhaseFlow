#!/usr/bin/env python3
"""Figure and LaTeX tables for the closed-cavity studies written by runStudy.

Usage: plot_cavity_results.py STUDY_DIR FIGURE.pdf TABLES.tex
"""
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Fixed categorical order, one colour per update in both panels.
UPDATES = {
    "backward": ("backward", "#2a78d6", "-", None),
    "blockTheta1": (r"monolithic $\theta=1$", "#eb6834", "none", "o"),
    "blockTheta0.55": (r"monolithic $\theta=0.55$", "#1baf7a", "-", None),
    "blockDirect": (r"monolithic $\theta=1/2$, direct", "#eda100", "-", None),
    "blockIterative": (r"monolithic $\theta=1/2$, iterative", "#e87ba4", "--",
                       None),
    "segregatedTheta": (r"segregated $\theta=1/2$", "#008300", ":", None),
}


def summary(case):
    return json.loads(
        (case / "postProcessing/cavityVerification/summary.json").read_text()
    )


def history(case):
    return np.loadtxt(case / "postProcessing/cavityVerification/history.csv",
                      delimiter=",", skiprows=1)


def modulus(s):
    xi, theta = s["xi"], s["setup"]["theta"]
    if s["setup"]["method"] == "fullyImplicit":
        return (1.0 + xi * xi) ** -0.5
    return math.sqrt((1 + (1 - theta) ** 2 * xi * xi) / (1 + theta**2 * xi * xi))


def energy_swing(case, s):
    """Relative peak-to-peak variation of E_ac from a fit at 2 omega_fit.

    A fit at the known frequency avoids the sampling bias of max - min.
    """
    h = history(case)
    t, energy = h[:, 0], h[:, 6]
    w = 2.0 * s["rates"]["omegaFit"]
    basis = np.column_stack((np.ones_like(t), np.cos(w * t), np.sin(w * t)))
    (c0, c1, c2), *_ = np.linalg.lstsq(basis, energy, rcond=None)
    return 2.0 * math.hypot(c1, c2) / c0


def sci(x, digits=2):
    if x is None:
        return "--"
    if x == 0:
        return "$0$"
    exponent = math.floor(math.log10(abs(x)))
    mantissa = x / 10**exponent
    if round(mantissa, digits - 1) >= 10:
        mantissa, exponent = mantissa / 10, exponent + 1
    return f"${mantissa:.{digits - 1}f}\\times10^{{{exponent}}}$"


def order(coarse, fine):
    return f"{math.log(coarse / fine) / math.log(2.0):.2f}"


def figure(study, path):
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
        "axes.linewidth": 0.6, "axes.edgecolor": "#555555",
        "xtick.color": "#555555", "ytick.color": "#555555",
        "axes.labelcolor": "#222222", "legend.frameon": False,
    })
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.3, 3.1))
    for name in ("backward", "blockTheta1", "blockTheta0.55", "blockDirect"):
        label, colour, style, marker = UPDATES[name]
        s, h = summary(study / name), history(study / name)
        periods = h[:, 0] * s["setup"]["f"]
        amp = np.sqrt(h[:, 5] / h[0, 5])
        every = s["setup"]["stepsPerPeriod"] * 4
        ax_a.semilogy(periods, amp, color=colour, lw=1.6 if marker is None else 0,
                      ls=style, marker=marker, ms=4, markevery=every,
                      mfc="none", mew=1.0, label=label)
        if name in ("backward", "blockTheta0.55"):
            n = np.arange(len(h))
            ax_a.semilogy(periods, modulus(s) ** n, color="#222222", lw=0.7,
                          ls=(0, (3, 2)))
    ax_a.plot([], [], color="#222222", lw=0.7, ls=(0, (3, 2)),
              label=r"$|g|^n$, analysis")
    ax_a.set_ylim(1e-18, 3)
    ax_a.set_xlim(0, 100)
    ax_a.set_xlabel(r"$t/T$")
    ax_a.set_ylabel(r"$(\mathcal{E}_h^n/\mathcal{E}_h^0)^{1/2}$")
    ax_a.set_title("(a) modal amplitude", loc="left", fontsize=9)

    for name in ("blockDirect", "blockIterative", "segregatedTheta"):
        label, colour, style, _ = UPDATES[name]
        s, h = summary(study / name), history(study / name)
        periods = h[:, 0] * s["setup"]["f"]
        drift = np.abs(h[:, 5] / h[0, 5] - 1.0)
        drift[0] = np.nan
        ax_b.semilogy(periods, drift, color=colour,
                      lw=1.4 if style == "-" else 1.0, ls=style, label=label)
    ax_b.set_xlim(0, 100)
    ax_b.set_ylim(1e-17, 1e-11)
    ax_b.set_xlabel(r"$t/T$")
    ax_b.set_ylabel(r"$|\mathcal{E}_h^n/\mathcal{E}_h^0-1|$")
    ax_b.set_title(r"(b) discrete energy, $\theta=1/2$", loc="left",
                   fontsize=9)
    for ax in (ax_a, ax_b):
        ax.grid(True, which="major", color="#e6e6e6", lw=0.5)
        ax.set_axisbelow(True)
    handles, labels = {}, []
    for ax in (ax_a, ax_b):
        for handle, label in zip(*ax.get_legend_handles_labels()):
            if label not in handles:
                handles[label] = handle
                labels.append(label)
    fig.legend([handles[l] for l in labels], labels, loc="lower center",
               ncol=3, fontsize=7.5, handlelength=2.6, columnspacing=1.4)
    fig.tight_layout(w_pad=1.5, rect=(0, 0.17, 1, 1))
    fig.savefig(path)


def tables(root, path):
    lines = []
    add = lines.append

    add("% Study A")
    rows = [
        ("backward", "backward", "--"),
        ("blockTheta1", r"monolithic, $\theta=1$", "direct"),
        ("blockTheta0.55", r"monolithic, $\theta=0.55$", "direct"),
        ("blockDirect", r"monolithic, $\theta=1/2$", "direct"),
        ("blockIterative", r"monolithic, $\theta=1/2$", "iterative"),
        ("segregatedTheta", r"segregated, $\theta=1/2$", "--"),
    ]
    for name, label, solver in rows:
        s = summary(root / "A" / name)
        gamma_fit = s["rates"]["gammaFit"]
        gamma_fit = (r"$<10^{-11}$" if abs(gamma_fit) < 1e-9
                     else f"${gamma_fit:.3f}$")
        gamma_disc = s["rates"]["gammaDiscrete"]
        gamma_disc = "$0$" if gamma_disc == 0 else f"${gamma_disc:.3f}$"
        drift = (sci(s["energy"]["discreteWaveEnergyMaxRelativeDrift"])
                 if s["setup"]["theta"] == 0.5
                 and s["setup"]["method"] != "fullyImplicit" else "--")
        add(f"    {label} & {solver} & "
            f"{sci(s['implementation']['maxModalErrorVsDiscrete'])} & "
            f"{gamma_fit} & {gamma_disc} & {drift} \\\\")

    add("% Study B, trapezoidal")
    trap = sorted((summary(p) for p in (root / "B").glob("trapezoidal_*")),
                  key=lambda s: s["setup"]["stepsPerPeriod"])
    for i, s in enumerate(trap):
        f = abs(s["rates"]["relativeFrequencyErrorVsSemiDiscrete"])
        e = s["E2"]["vsSemiDiscrete"]
        prev = trap[i - 1] if i else None
        qf = order(abs(prev["rates"]["relativeFrequencyErrorVsSemiDiscrete"]),
                   f) if prev else "--"
        qe = order(prev["E2"]["vsSemiDiscrete"], e) if prev else "--"
        add(f"    ${s['setup']['stepsPerPeriod']}$ & ${s['xi']:.4f}$ & "
            f"${s['courant']:g}$ & {sci(f)} & {qf} & {sci(e)} & {qe} \\\\")

    add("% Study B, backward")
    back = sorted((summary(p) for p in (root / "B").glob("backward_*")),
                  key=lambda s: s["setup"]["stepsPerPeriod"])
    for i, s in enumerate(back):
        f = abs(s["rates"]["relativeFrequencyErrorVsSemiDiscrete"])
        d = s["rates"]["gammaFit"] / s["omega_h"]
        e = s["E2"]["vsSemiDiscrete"]
        amp = s["endOfRun"]["amplitudeFactorFit"]
        prev = back[i - 1] if i else None
        q = [order(abs(prev["rates"]["relativeFrequencyErrorVsSemiDiscrete"]),
                   f),
             order(prev["rates"]["gammaFit"] / prev["omega_h"], d),
             order(prev["E2"]["vsSemiDiscrete"], e)] if prev else ["--"] * 3
        add(f"    ${s['setup']['stepsPerPeriod']}$ & {sci(f)} & {q[0]} & "
            f"{sci(d)} & {q[1]} & {sci(amp, 3)} & {sci(e)} & {q[2]} \\\\")

    add("% Study C")
    cases = sorted((root / "C").glob("trapezoidal_*"),
                   key=lambda p: summary(p)["setup"]["cellsPerWavelength"])
    spat = [summary(p) for p in cases]
    for i, s in enumerate(spat):
        kh = s["setup"]["k"] * s["setup"]["h"]
        rec = abs(s["rates"]["relativeSpatialFrequencyError"])
        pred = abs(s["rates"]["relativeSpatialFrequencyErrorPredicted"])
        ep, eu = s["E2"]["vsExact"], s["E2"]["velocityVsExact"]
        swing = energy_swing(cases[i], s)
        prev = spat[i - 1] if i else None
        q = [order(abs(prev["rates"]["relativeSpatialFrequencyError"]), rec),
             order(prev["E2"]["vsExact"], ep),
             order(prev["E2"]["velocityVsExact"], eu)] if prev else ["--"] * 3
        add(f"    ${s['setup']['cellsPerWavelength']}$ & ${kh:.4f}$ & "
            f"{sci(rec, 4)} & {sci(pred, 4)} & {q[0]} & {sci(ep)} & {q[1]} & "
            f"{sci(eu)} & {q[2]} & ${swing / (0.5 * kh) ** 2:.3f}$ \\\\")

    Path(path).write_text("\n".join(lines) + "\n")


def main():
    root = Path(sys.argv[1])
    figure(root / "A", sys.argv[2])
    tables(root, sys.argv[3])


if __name__ == "__main__":
    main()
