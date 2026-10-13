#!/usr/bin/env python3
"""Figure and LaTeX table rows for the amplification-matrix study C.

Usage: plot_spectral_results.py STUDY_C_DIR FIGURE.pdf TABLES.tex
"""
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402

UPDATES = {
    "monolithic": ("monolithic", "#0b2e6b", "o"),
    "segregated": ("segregated", "#c4402f", "s"),
    "segregatedImplicit": ("segregated, implicit damping", "#2a9d6f", "^"),
}
# Horizontal offsets so that coinciding markers at rho = 1 stay visible.
OFFSET = {"monolithic": 0.92, "segregated": 1.0, "segregatedImplicit": 1.08}
COURANTS = ("0.5", "2", "8")
SIGMAS = ("1.26e5", "2.53e5", "5.05e5")
FILL = {"0.5": "none", "2": "half", "8": "full"}


def spectrum(case):
    return json.loads((case / "spectrum.json").read_text())


def sci(x, digits=2):
    if x <= 0:
        return "$0$"
    exponent = math.floor(math.log10(abs(x)))
    mantissa = x / 10**exponent
    if round(mantissa, digits - 1) >= 10:
        mantissa, exponent = mantissa / 10, exponent + 1
    return f"${mantissa:.{digits - 1}f}\\times10^{{{exponent}}}$"


def figure(root, path):
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
        "axes.linewidth": 0.6, "axes.edgecolor": "#555555",
        "xtick.color": "#555555", "ytick.color": "#555555",
        "axes.labelcolor": "#222222", "legend.frameon": False,
    })
    fig, ax = plt.subplots(figsize=(5.2, 3.3))
    x = np.geomspace(0.15, 25, 400)
    ax.plot(x, np.maximum(1, np.abs(1 - x)), color="#999999", lw=0.9,
            ls="--", zorder=1)
    for key, (label, colour, marker) in UPDATES.items():
        for co in COURANTS:
            for sigma in SIGMAS:
                r = spectrum(root / f"{key}_Co{co}_sigma{sigma}")
                fill = FILL[co]
                ax.plot(r["outerCellSigmaDt"] * OFFSET[key], r["spectralRadius"],
                        ls="none", marker=marker, ms=6, color=colour,
                        mfc=("none" if fill == "none" else colour),
                        fillstyle=("full" if fill != "half" else "left"),
                        mew=1.0, zorder=3)
    scan = sorted(root.glob("scan_sigmaDt*"))
    ax.plot([spectrum(c)["outerCellSigmaDt"] for c in scan],
            [spectrum(c)["spectralRadius"] for c in scan], ls="none",
            marker="x", ms=5, color=UPDATES["segregated"][1], mew=1.0,
            zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.15, 25)
    ax.set_ylim(0.8, 20)
    ax.set_xlabel(r"$\sigma\Delta t$ in the outermost PML cell")
    ax.set_ylabel(r"spectral radius $\rho(\mathsf{G})$")
    ax.grid(True, color="#e6e6e6", lw=0.5)
    ax.set_axisbelow(True)
    handles = [Line2D([], [], ls="none", marker=m, color=c, mfc=c, ms=6,
                      label=l) for l, c, m in UPDATES.values()]
    handles.append(Line2D([], [], ls="none", marker="x", color=UPDATES["segregated"][1],
                          ms=5, label=r"segregated, $\mathrm{Co}=2$ scan"))
    handles.append(Line2D([], [], color="#999999", ls="--",
                          label=r"$\max(1,|1-\sigma\Delta t|)$"))
    handles += [Line2D([], [], ls="none", marker="o", color="#444444", ms=6,
                       mfc=("none" if FILL[c] == "none" else "#444444"),
                       fillstyle=("full" if FILL[c] != "half" else "left"),
                       label=fr"$\mathrm{{Co}}={c}$") for c in COURANTS]
    ax.legend(handles=handles, fontsize=7, loc="upper left", ncol=1)
    fig.tight_layout()
    fig.savefig(path)


def tables(root, path):
    lines = ["% Main"]
    for co in COURANTS:
        for sigma in SIGMAS:
            m = spectrum(root / f"monolithic_Co{co}_sigma{sigma}")
            s = spectrum(root / f"segregated_Co{co}_sigma{sigma}")
            i = spectrum(root / f"segregatedImplicit_Co{co}_sigma{sigma}")
            cells = [f"${co}$", f"${float(sigma) / 1e5:.2f}$",
                     f"${m['outerCellSigmaDt']:.3g}$",
                     sci(m["spectralRadius"] - 1, 1),
                     f"${m['powerNorm']['last'] / 1e3:.1f}\\times10^{{3}}$",
                     f"${s['spectralRadius']:.5f}$",
                     f"${s['explicitDampingRadius']:.5f}$",
                     f"${i['spectralRadius']:.4f}$"]
            lines.append("    " + " & ".join(cells) + " \\\\")
    lines.append("% Scan")
    for case in sorted(root.glob("scan_sigmaDt*")):
        r = spectrum(case)
        lines.append(f"    ${r['outerCellSigmaDt']:.2f}$ & "
                     f"${r['spectralRadius']:.5f}$ & ${r['nUnstable']}$ \\\\")
    Path(path).write_text("\n".join(lines) + "\n")


def main():
    root = Path(sys.argv[1])
    figure(root, sys.argv[2])
    tables(root, sys.argv[3])


if __name__ == "__main__":
    main()
