#!/usr/bin/env python3
"""Figure and LaTeX table rows for the radiation-average studies.

Usage: plot_radiation_results.py CAVITY_DIR TWO_LAYER_DIR FIGURE.pdf TABLES.tex

CAVITY_DIR and TWO_LAYER_DIR are the studies/radiation directories of
closedCavity1D (studies A and B) and twoLayerCavity1D (study C).
"""
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def summary(case):
    return json.loads((case / "postProcessing/radiationVerification/"
                       "summary.json").read_text())


def sci(x, digits=2):
    exponent = math.floor(math.log10(abs(x)))
    mantissa = x / 10**exponent
    if round(mantissa, digits - 1) >= 10:
        mantissa, exponent = mantissa / 10, exponent + 1
    return f"${mantissa:.{digits - 1}f}\\times10^{{{exponent}}}$"


def order(coarse, fine):
    return f"{math.log(coarse / fine) / math.log(2.0):.2f}"


def whole(m0):
    return abs(2 * m0 - round(2 * m0)) < 1e-9


def figure(items, path):
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
        "axes.linewidth": 0.6, "axes.edgecolor": "#555555",
        "xtick.color": "#555555", "ytick.color": "#555555",
        "axes.labelcolor": "#222222", "legend.frameon": False,
    })
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    m = np.array([s["setup"]["radiationPeriods"] for s in items])
    err = np.array([s["errors"]["mean"]["pr"] for s in items])
    bound = np.array([s["biasBoundDiscrete"] for s in items])
    w = np.array([whole(v) for v in m])
    m_line = np.geomspace(0.4, 30, 200)
    ax.plot(m_line, 1 / (2 * np.pi * m_line), color="#999999", lw=0.8,
            ls="--", label=r"$1/(2\pi M_0)$")
    ax.plot(m[~w], bound[~w], ls="none", marker="x", ms=5, color="#222222",
            label="bias bound")
    ax.plot(m[~w], err[~w], ls="none", marker="s", ms=4, mfc="none",
            mew=1.0, color="#2a78d6", label=r"$2M_0$ not an integer")
    ax.plot(m[w], err[w], ls="none", marker="o", ms=4, color="#0b2e6b",
            label=r"$2M_0$ an integer")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"window length $M_0$ (periods)")
    ax.set_ylabel(r"error of $p_r$ against the exact mean")
    ax.grid(True, color="#e6e6e6", lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(path)


def tables(cavity, two_layer, path):
    lines = []
    add = lines.append

    add("% Window")
    items = sorted((summary(c) for c in (cavity / "A").iterdir()),
                   key=lambda s: s["setup"]["radiationPeriods"])
    whole_items = [s for s in items if whole(s["setup"]["radiationPeriods"])]
    other_items = [s for s in items if not whole(s["setup"]["radiationPeriods"])]
    for a, b in zip(whole_items, other_items):
        cells = []
        for s in (a, b):
            cells += [f"${s['setup']['radiationPeriods']:g}$",
                      f"${s['windowSamples']}$",
                      sci(s["biasBoundDiscrete"]),
                      sci(s["errors"]["mean"]["pr"]),
                      sci(s["errors"]["mean"]["Rxx"])]
        add("    " + " & ".join(cells) + " \\\\")
    discrete = max(max(s["errors"]["discrete"].values()) for s in items)
    add(f"% max error against the discrete reference: {discrete:.2e}")

    add("% Cavity refinement")
    items = sorted((summary(c) for c in (cavity / "B").iterdir()),
                   key=lambda s: s["setup"]["cellsPerWavelength"])
    for i, s in enumerate(items):
        e, prev = s["errors"], items[i - 1]["errors"] if i else None
        q = [order(prev["mean"][k], e["mean"][k]) if i else "--"
             for k in ("pr", "Rxx")]
        add(f"    ${s['setup']['cellsPerWavelength']}$ & "
            f"${s['setup']['stepsPerPeriod']}$ & "
            f"{sci(max(e['discrete'].values()), 1)} & "
            f"{sci(e['mean']['pr'])} & {q[0]} & "
            f"{sci(e['mean']['Rxx'])} & {q[1]} \\\\")

    add("% Two-layer refinement (pr)")
    groups = ("gasInterior", "gasInterfaceCell", "liquidInterfaceCell",
              "liquidInterior")
    for mode in ("lower", "upper"):
        items = sorted((summary(c) for c in (two_layer / "C").glob(f"{mode}_*")),
                       key=lambda s: s["setup"]["cellsPerUnit"])
        for i, s in enumerate(items):
            e, prev = s["errors"]["mean"], items[i - 1]["errors"]["mean"] if i else None
            disc = max(v for g in s["errors"]["discrete"].values()
                       for v in g.values())
            cells = [mode, f"${s['setup']['nx']}$",
                     f"${s['setup']['stepsPerPeriod']}$", sci(disc, 1)]
            for g in groups:
                cells += [sci(e[g]["pr"]),
                          order(prev[g]["pr"], e[g]["pr"]) if i else "--"]
            add("    " + " & ".join(cells) + " \\\\")

    Path(path).write_text("\n".join(lines) + "\n")


def main():
    cavity, two_layer = Path(sys.argv[1]), Path(sys.argv[2])
    items = sorted((summary(c) for c in (cavity / "A").iterdir()),
                   key=lambda s: s["setup"]["radiationPeriods"])
    figure(items, sys.argv[3])
    tables(cavity, two_layer, sys.argv[4])


if __name__ == "__main__":
    main()
