#!/usr/bin/env python3
"""Figure and LaTeX table rows for the two-layer cavity studies.

Usage: plot_two_layer_results.py STUDY_DIR FIGURE.pdf TABLES.tex
"""
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Fixed order: one colour per interface position, one marker per mode.
INTERFACES = [
    (0.0, "on a face", "#2a78d6", "-"),
    (0.5, r"cut, $s=1/2$", "#1baf7a", "--"),
    (0.25, r"cut, $s=1/4$", "#eb6834", ":"),
]
MODES = [("lower", "o", "lower mode"), ("upper", "s", "upper mode")]


def summary(case):
    return json.loads((case / "postProcessing/twoLayerVerification/"
                       "summary.json").read_text())


def probe(case, field):
    rows = [line.replace("(", " ").replace(")", " ").split()
            for line in (case / f"postProcessing/cellProbes/0/{field}")
            .read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    data = np.array(rows, dtype=float)
    return data[:, 1:] if field == "p" else data[:, 1::3]


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


def refinement(root, letter, offset, mode):
    cases = [c for c in (root / letter).iterdir()
             if summary(c)["setup"]["interfaceOffset"] == offset
             and summary(c)["setup"]["mode"] == mode]
    return sorted((summary(c) for c in cases),
                  key=lambda s: s["setup"]["cellsPerUnit"])


def figure(root, path):
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
        "axes.linewidth": 0.6, "axes.edgecolor": "#555555",
        "xtick.color": "#555555", "ytick.color": "#555555",
        "axes.labelcolor": "#222222", "legend.frameon": False,
    })
    fig, ax = plt.subplots(figsize=(4.6, 3.3))
    for offset, label, colour, style in INTERFACES:
        letter = "B" if offset == 0 else "C"
        for mode, marker, _ in MODES:
            rows = refinement(root, letter, offset, mode)
            n = [32 * s["setup"]["cellsPerUnit"] for s in rows]
            err = [abs(s["rates"]["relativeFrequencyErrorFitVsExact"])
                   for s in rows]
            ax.loglog(n, err, color=colour, ls=style, lw=1.3, marker=marker,
                      ms=5, mfc="white" if mode == "upper" else colour,
                      mew=1.1)
    for slope, anchor in ((1, 2.5e-3), (2, 1.2e-5)):
        x = np.array([128, 512])
        ax.loglog(x, anchor * (x / 512.0) ** (-slope), color="#888888",
                  lw=0.8)
        ax.text(560, anchor * 0.95, f"order {slope}", fontsize=7.5,
                color="#555555", va="center")
    handles = [plt.Line2D([], [], color=c, ls=s, lw=1.3, label=l)
               for _, l, c, s in INTERFACES]
    handles += [plt.Line2D([], [], color="#555555", ls="none", marker=m,
                           mfc="white" if mode == "upper" else "#555555",
                           label=l) for mode, m, l in MODES]
    ax.legend(handles=handles, fontsize=7.5, loc="lower left", ncol=2)
    ax.set_xlabel(r"gas cells per wavelength $N_{\lambda,g}$")
    ax.set_ylabel(r"$|\omega_h/\omega-1|$")
    ax.set_xlim(25, 900)
    ax.set_xticks([32, 64, 128, 256, 512])
    ax.set_xticklabels(["32", "64", "128", "256", "512"])
    ax.minorticks_off()
    ax.grid(True, which="major", color="#e6e6e6", lw=0.5)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(path)


def tables(root, path):
    lines = []
    add = lines.append

    add("% Study A")
    labels = {"backward": ("backward", "--"),
              "segregated": (r"segregated, $\theta=1/2$", "--"),
              "blockDirect": (r"monolithic, $\theta=1/2$", "direct"),
              "blockIterative": (r"monolithic, $\theta=1/2$", "iterative")}
    for mode, _, _ in MODES:
        for key, (label, solver) in labels.items():
            s = summary(root / "A" / f"{mode}_{key}")
            g = s["rates"]["gammaFit"]
            g = r"$<10^{-10}$" if abs(g) < 1e-9 else f"${g:.2f}$"
            gd = s["rates"]["gammaDiscrete"]
            gd = r"$<10^{-10}$" if abs(gd) < 1e-9 else f"${gd:.2f}$"
            drift = (sci(s["energy"]["discreteWaveEnergyMaxRelativeDrift"])
                     if s["setup"]["method"] != "fullyImplicit" else "--")
            add(f"    {mode} & {label} & {solver} & "
                f"{sci(s['implementation']['maxModalErrorVsDiscrete'])} & "
                f"{g} & {gd} & {drift} \\\\")

    for letter, offsets in (("B", (0.0,)), ("C", (0.25, 0.5))):
        add(f"% Study {letter}")
        for offset in offsets:
            for mode, _, _ in MODES:
                rows = refinement(root, letter, offset, mode)
                for i, s in enumerate(rows):
                    f = abs(s["rates"]["relativeFrequencyErrorFitVsExact"])
                    ep, eu = s["E2"]["pressureVsExact"], s["E2"]["velocityVsExact"]
                    if i:
                        p = rows[i - 1]
                        q = [order(abs(p["rates"]
                                       ["relativeFrequencyErrorFitVsExact"]), f),
                             order(p["E2"]["pressureVsExact"], ep),
                             order(p["E2"]["velocityVsExact"], eu)]
                    else:
                        q = ["--"] * 3
                    lead = (f"$1/{round(1 / offset)}$ & " if letter == "C"
                            else "")
                    add(f"    {lead}{mode} & {32 * s['setup']['cellsPerUnit']} & "
                        f"{sci(f)} & {q[0]} & {sci(ep)} & {q[1]} & "
                        f"{sci(eu)} & {q[2]} \\\\")

    add("% Study D")
    for offset, label in ((0.0, "face"), (0.5, "$s=1/2$")):
        for mode, _, _ in MODES:
            tag = "s0" if offset == 0 else "s0.5"
            serial = root / "D" / f"{tag}_{mode}_np1"
            parallel = root / "D" / f"{tag}_{mode}_np2"
            s2 = summary(parallel)
            proc = s2["processorFace"]
            dp = float(np.abs(probe(serial, "p") - probe(parallel, "p")).max())
            du = float(np.abs(probe(serial, "U") - probe(parallel, "U")).max())
            add(f"    {label} & {mode} & ${proc['invRhofRank0']:.5g}$ & "
                f"${proc['invRhofRank1']:.5g}$ & "
                f"{sci(s2['energy']['discreteWaveEnergyMaxRelativeDrift'])} & "
                f"{sci(dp)} & {sci(du)} \\\\")

    Path(path).write_text("\n".join(lines) + "\n")


def main():
    root = Path(sys.argv[1])
    figure(root, sys.argv[2])
    tables(root, sys.argv[3])


if __name__ == "__main__":
    main()
