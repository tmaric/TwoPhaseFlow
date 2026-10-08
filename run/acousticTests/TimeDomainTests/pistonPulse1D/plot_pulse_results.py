#!/usr/bin/env python3
"""Figure and LaTeX table rows for the piston-pulse studies.

Usage: plot_pulse_results.py STUDY_DIR FIGURE.pdf TABLES.tex
"""
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pulseReference import LayeredPulseSolution, PistonPulse  # noqa: E402
from summarize_pulse_study import temporal_errors  # noqa: E402


def summary(case):
    return json.loads((case / "postProcessing/pulseVerification/"
                       "summary.json").read_text())


def sci(x, digits=2):
    if x is None:
        return "--"
    exponent = math.floor(math.log10(abs(x)))
    mantissa = x / 10**exponent
    if round(mantissa, digits - 1) >= 10:
        mantissa, exponent = mantissa / 10, exponent + 1
    return f"${mantissa:.{digits - 1}f}\\times10^{{{exponent}}}$"


def order(coarse, fine):
    return f"{math.log(coarse / fine) / math.log(2.0):.2f}"


def figure(case, path):
    s = json.loads((case / "caseSetup.json").read_text())
    trace = np.loadtxt(case / "postProcessing/pulseVerification/traces.csv",
                       delimiter=",", skiprows=1)
    names = (case / "postProcessing/pulseVerification/traces.csv"
             ).read_text().splitlines()[0].split(",")
    pulse = PistonPulse(s["f0"], s["nCycles"], s["u0"])
    exact = LayeredPulseSolution(pulse, s["rho1"], s["c1"], s["rho2"],
                                 s["c2"], s["xInterface"])
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
        "axes.linewidth": 0.6, "axes.edgecolor": "#555555",
        "xtick.color": "#555555", "ytick.color": "#555555",
        "axes.labelcolor": "#222222", "legend.frameon": False,
    })
    fig, axes = plt.subplots(2, 1, figsize=(6.3, 3.9), sharex=True)
    t_fine = np.linspace(0.0, s["endTime"], 6000)
    for ax, column, label in ((axes[0], 1, "(a) gas, $x=0.045\\,$m"),
                              (axes[1], 7, "(b) liquid, $x=0.145\\,$m")):
        x = float(names[column].split("_x")[1])
        p_line, _ = exact.fields([x], t_fine)
        ax.plot(t_fine * 1e3, p_line[:, 0], color="#222222", lw=0.9,
                label="exact")
        ax.plot(trace[:, 0] * 1e3, trace[:, column], ls="none", marker="o",
                ms=2.2, color="#2a78d6", label="computed")
        ax.set_ylabel(r"$p$ (Pa)")
        ax.set_title(label, loc="left", fontsize=9)
        ax.grid(True, color="#e6e6e6", lw=0.5)
        ax.set_axisbelow(True)
    # The liquid probe is quiet between its first two passages.
    axes[1].legend(fontsize=7.5, loc="upper left", ncol=2,
                   bbox_to_anchor=(0.31, 1.0))
    axes[1].set_xlabel(r"$t$ (ms)")
    axes[1].set_xlim(0, s["endTime"] * 1e3)
    fig.tight_layout()
    fig.savefig(path)


def tables(root, path):
    lines = []
    add = lines.append

    add("% Temporal")
    trap = {summary(c)["setup"]["stepsPerPeriod"]: (summary(c), c)
            for c in (root / "A").glob("trapezoidal_*")}
    back = {summary(c)["setup"]["stepsPerPeriod"]: (summary(c), c)
            for c in (root / "A").glob("backward_*")}
    t_err = temporal_errors(root, [c for _, c in trap.values()])
    b_err = temporal_errors(root, [c for _, c in back.values()])
    previous = {}
    for nt in sorted(set(trap) | set(back)):
        cells = [f"${nt}$", f"${2 * math.pi / nt:.4f}$",
                 f"${137.2 / nt:.3g}$"]
        if nt in trap:
            s, c = trap[nt]
            e = t_err[c.name]
            q = order(previous["t"], e) if "t" in previous else "--"
            cells += [sci(s["E2"]["pressure"]), sci(e), q]
            previous["t"] = e
        else:
            cells += ["", "", ""]
        s, c = back[nt]
        e = b_err[c.name]
        q = order(previous["b"], e) if "b" in previous else "--"
        cells += [sci(e), q]
        previous["b"] = e
        add("    " + " & ".join(cells) + " \\\\")

    add("% Spatial")
    for letter, label in (("B", "gas first"), ("D", "liquid first")):
        items = sorted((summary(c) for c in (root / letter).iterdir()),
                       key=lambda s: s["setup"]["nx"])
        for i, s in enumerate(items):
            e = s["E2"]
            p = items[i - 1]["E2"] if i else None
            q = [order(p[k], e[k]) if i else "--"
                 for k in ("pressure", "velocity", "pressureMedium2")]
            spur = s["pml"]["spuriousRelativePressure"]
            add(f"    {label} & ${s['setup']['gasCellsPerWavelength']:.0f}$ & "
                f"${s['setup']['stepsPerPeriod']}$ & {sci(e['pressure'])} & "
                f"{q[0]} & {sci(e['velocity'])} & {q[1]} & "
                f"{sci(e['pressureMedium2'])} & {q[2]} & "
                f"{sci(spur) if spur else '--'} \\\\")

    add("% PML")
    items = sorted((summary(c) for c in (root / "C").iterdir()),
                   key=lambda s: (s["setup"]["pmlThickness"],
                                  s["setup"]["sigmaMax"]))
    for s in items:
        st = s["setup"]
        factor = math.exp(-2 * st["sigmaMax"] * st["pmlThickness"]
                          / ((st["pmlOrder"] + 1) * st["c2"]))
        spur = s["pml"]["spuriousRelativePressure"]
        # The ratio is meaningful only while the layer, not the mesh-dependent
        # floor (about 3e-5 here), dominates the measured pressure.
        ratio = f"${spur / factor:.2f}$" if factor > 4e-5 else "--"
        add(f"    ${st['pmlThickness']:g}$ & {sci(st['sigmaMax'], 1)} & "
            f"{sci(factor)} & {sci(spur)} & {ratio} \\\\")

    Path(path).write_text("\n".join(lines) + "\n")


def main():
    root = Path(sys.argv[1])
    figure(root / "B" / "NX2800", sys.argv[2])
    tables(root, sys.argv[3])


if __name__ == "__main__":
    main()
