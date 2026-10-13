#!/usr/bin/env python3
"""Figure and LaTeX table rows for the multidirectional PML studies.

Usage: plot_pml_results.py STUDY_DIR FIGURE.pdf TABLES.tex
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summarize_pml_study import pressure, T_OBS  # noqa: E402

COURANTS = ("0.5", "2", "8")
SIGMAS = ("1.26e5", "2.53e5", "5.05e5")
COLOURS = {"0.5": "#7fb2ea", "2": "#2a78d6", "8": "#0b2e6b"}
STYLES = {"1.26e5": ":", "2.53e5": "--", "5.05e5": "-"}
CAP = 1e3


def summary(case):
    return json.loads((case / "postProcessing/pmlVerification/"
                       "summary.json").read_text())


def history(case, name):
    path = case / "postProcessing/pmlVerification" / name
    return np.loadtxt(path, delimiter=",", skiprows=1)


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


def figure(root, path):
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
        "axes.linewidth": 0.6, "axes.edgecolor": "#555555",
        "xtick.color": "#555555", "ytick.color": "#555555",
        "axes.labelcolor": "#222222", "legend.frameon": False,
    })
    fig, axes = plt.subplots(2, 2, figsize=(6.3, 4.6), sharex=True,
                             sharey="row")
    titles = {"blockCoupled": "monolithic", "thetaNewmark": "segregated"}
    for col, method in enumerate(("blockCoupled", "thetaNewmark")):
        for co in COURANTS:
            for sigma in SIGMAS:
                case = root / "B" / f"{method}_Co{co}_sigma{sigma}"
                energy = history(case, "energy.csv")
                maxima = history(case, "maxima.csv")
                crashed = summary(case)["lastSampleTime"] < 0.049
                for row, (t, y) in enumerate((
                        (energy[:, 0], energy[:, 1] / energy[0, 1]),
                        (maxima[:, 0], maxima[:, 1]))):
                    ax = axes[row, col]
                    keep = np.isfinite(y) & (y < CAP) & (t > 0)
                    ax.plot(t[keep] * 1e3, y[keep], color=COLOURS[co],
                            ls=STYLES[sigma], lw=1.0)
                    if crashed:
                        last = np.flatnonzero(keep)[-1]
                        ax.plot(t[last] * 1e3, y[last], marker="x", ms=5,
                                mew=1.2, color=COLOURS[co])
        axes[0, col].set_title(f"({'ab'[col]}) {titles[method]}", loc="left",
                               fontsize=9)
        axes[1, col].set_title(f"({'cd'[col]}) {titles[method]}", loc="left",
                               fontsize=9)
        axes[1, col].set_xlabel(r"$t$ (ms)")
    axes[0, 0].set_ylabel(r"$E_{\mathrm{ac}}(t)/E_{\mathrm{ac}}(0)$")
    axes[1, 0].set_ylabel(r"$\max|p|$ (Pa), whole domain")
    for ax in axes.flat:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(0.08, 50)
        ax.grid(True, color="#e6e6e6", lw=0.5)
        ax.set_axisbelow(True)
    axes[0, 0].set_ylim(1e-12, CAP)
    axes[1, 0].set_ylim(1e-9, CAP)
    handles = [Line2D([], [], color=COLOURS[c], lw=1.4,
                      label=fr"$\mathrm{{Co}}={c}$") for c in COURANTS]
    handles += [Line2D([], [], color="#444444", ls=STYLES[s], lw=1.0,
                       label=fr"$\sigma_{{\max}}={float(s) / 1e5:.2f}\times10^5$")
                for s in SIGMAS]
    handles.append(Line2D([], [], color="#444444", ls="none", marker="x",
                          label=r"exceeds $10^3$"))
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7.5,
               bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(path)


def tables(root, path):
    lines = []
    add = lines.append

    add("% Accuracy")
    previous = None
    for h in ("2e-3", "1e-3", "5e-4"):
        pml, ref = root / "A" / f"pml_h{h}", root / "A" / f"enlarged_h{h}"
        s, r = summary(pml), summary(ref)
        t, p = pressure(pml)
        _, p_ref = pressure(ref)
        obs = t <= T_OBS * (1 + 1e-9)
        trunc = math.sqrt(np.sum((p[obs] - p_ref[obs]) ** 2)
                          / np.sum(p_ref[obs] ** 2))
        e = (r["E2exact"]["pressure"], r["E2exact"]["velocity"],
             s["E2exact"]["pressure"], trunc)
        q = [order(a, b) if previous else "--" for a, b in
             zip(previous or e, e)]
        add(f"    ${s['setup']['h'] * 1e3:g}$ & "
            f"${s['setup']['pulseWidth'] / s['setup']['h']:g}$ & "
            + " & ".join(f"{sci(v)} & {o}" for v, o in zip(e, q)) + " \\\\")
        previous = e

    add("% Long time")
    for co in COURANTS:
        for sigma in SIGMAS:
            cells = []
            for name in (f"blockCoupled_Co{co}_sigma{sigma}",
                         f"blockCoupled_L0.012_Co{co}_sigma{sigma}",
                         f"thetaNewmark_Co{co}_sigma{sigma}"):
                s = summary(root / "B" / name)
                if s["lastSampleTime"] < 0.049:
                    cells.append(f"stopped, ${s['lastSampleTime'] * 1e3:.1f}"
                                 "\\,\\mathrm{ms}$")
                else:
                    cells.append(sci(s["energy"]["final"]["Eac"]))
            st = s["setup"]
            add(f"    ${co}$ & ${float(sigma) / 1e5:.2f}$ & "
                f"${st['sigmaMax'] * st['deltaT']:.3g}$ & "
                + " & ".join(cells) + " \\\\")

    Path(path).write_text("\n".join(lines) + "\n")


def main():
    root = Path(sys.argv[1])
    figure(root, sys.argv[2])
    tables(root, sys.argv[3])


if __name__ == "__main__":
    main()
