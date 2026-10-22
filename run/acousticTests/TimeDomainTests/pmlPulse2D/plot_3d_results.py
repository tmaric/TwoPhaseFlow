#!/usr/bin/env python3
"""Figure and LaTeX table rows for the three-dimensional study D.

Usage: plot_3d_results.py STUDY_D_DIR FIGURE.pdf TABLES.tex

Panel (a): physical energy of the D2 runs. Panel (b): a Gaussian pulse on
the reduced mesh of D3 propagated with powers of G (monolithic and
segregated, Co 0.5, sigmaMax 5.05e5) to n = 2^24; the evolution is cached
in pulsePowers.json of each matrix directory. Table rows: D2 and D3.
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

MONO, SEG = "#0b2e6b", "#c4402f"


def sci(x, digits=2):
    if x is None or x <= 0:
        return "--"
    exponent = math.floor(math.log10(abs(x)))
    mantissa = x / 10**exponent
    if round(mantissa, digits - 1) >= 10:
        mantissa, exponent = mantissa / 10, exponent + 1
    return f"${mantissa:.{digits - 1}f}\\times10^{{{exponent}}}$"


def pulse_powers(case, kmax=24):
    cache = case / "pulsePowers.json"
    if cache.is_file():
        return json.loads(cache.read_text())
    data = np.load(case / "G.npz")
    g, names, n = data["G"], list(data["names"]), int(data["nCells"])
    s = json.loads((case / "setup.json").read_text())["setup"]
    h = s["h"]
    c = (np.arange(s["n"]) - s["n"] / 2 + 0.5) * h
    z, y, x = np.meshgrid(c, c, c, indexing="ij")
    # Pulse of width h near the centre of the physical region.
    p0 = np.exp(-((x - h / 2) ** 2 + (y - h / 2) ** 2 + z**2)
                / (2 * h**2)).ravel()
    y0 = np.zeros(len(g))
    y0[:n] = p0
    if "pPrev" in names:
        k = names.index("pPrev")
        y0[k * n:(k + 1) * n] = p0
    block = {nm: slice(i * n, (i + 1) * n) for i, nm in enumerate(names)}
    out = {"n": [], "p": [], "psi": [], "w": []}
    power = g.copy()
    for k in range(kmax + 1):
        state = power @ y0
        out["n"].append(2**k)
        out["p"].append(float(np.abs(state[block["p"]]).max()))
        out["psi"].append(float(np.abs(state[block["psi"]]).max()))
        out["w"].append(float(max(np.abs(state[block[w]]).max()
                                  for w in ("wx", "wy", "wz"))))
        power = power @ power
    out["deltaT"] = s["deltaT"]
    cache.write_text(json.dumps(out, indent=2) + "\n")
    return out


def figure(root, path):
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
        "axes.linewidth": 0.6, "axes.edgecolor": "#555555",
        "xtick.color": "#555555", "ytick.color": "#555555",
        "axes.labelcolor": "#222222", "legend.frameon": False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.9))
    ax = axes[0]
    for method, colour in (("blockCoupled", MONO), ("thetaNewmark", SEG)):
        for co, ls in (("0.5", "-"), ("2", "--")):
            case = root / "D2" / f"{method}_Co{co}"
            e = np.loadtxt(case / "postProcessing/pmlVerification/energy.csv",
                           delimiter=",", skiprows=1)
            y = e[:, 1] / e[0, 1]
            keep = (e[:, 0] > 0) & np.isfinite(y) & (y < 1e3)
            ax.plot(e[keep, 0] * 1e3, y[keep], color=colour, ls=ls, lw=1.0)
            if not keep[1:].all():
                last = np.flatnonzero(keep)[-1]
                ax.plot(e[last, 0] * 1e3, y[last], marker="x", ms=5,
                        mew=1.2, color=colour)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.08, 50)
    ax.set_ylim(1e-12, 1e3)
    ax.set_xlabel(r"$t$ (ms)")
    ax.set_ylabel(r"$E_{\mathrm{ac}}(t)/E_{\mathrm{ac}}(0)$")
    ax.set_title("(a) D2, $36^3$ cells", loc="left", fontsize=9)
    ax.legend(handles=[
        Line2D([], [], color=MONO, lw=1.2, label="monolithic"),
        Line2D([], [], color=SEG, lw=1.2, label="segregated"),
        Line2D([], [], color="#444444", ls="-", lw=1.0, label=r"$\mathrm{Co}=0.5$"),
        Line2D([], [], color="#444444", ls="--", lw=1.0, label=r"$\mathrm{Co}=2$"),
        Line2D([], [], color="#444444", ls="none", marker="x", label=r"exceeds $10^3$"),
    ], fontsize=7, loc="lower left")

    ax = axes[1]
    for name, colour in (("monolithic", MONO), ("segregated", SEG)):
        r = pulse_powers(root / "D3" / f"{name}_Co0.5_sigma5.05e5")
        n = np.array(r["n"])
        for key, ls in (("p", ":"), ("psi", "--"), ("w", "-")):
            ax.plot(n, np.maximum(r[key], 1e-30), color=colour, ls=ls, lw=1.0)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_ylim(1e-30, 1e2)
    ax.set_xlabel(r"time steps $n$")
    ax.set_ylabel("largest value, scaled state")
    ax.set_title(r"(b) D3, $10^3$ cells, $\mathrm{Co}=0.5$", loc="left",
                 fontsize=9)
    ax.legend(handles=[
        Line2D([], [], color=MONO, lw=1.2, label="monolithic"),
        Line2D([], [], color=SEG, lw=1.2, label="segregated"),
        Line2D([], [], color="#444444", ls=":", lw=1.0, label=r"$|p|$"),
        Line2D([], [], color="#444444", ls="--", lw=1.0, label=r"$|\psi|$"),
        Line2D([], [], color="#444444", ls="-", lw=1.0, label=r"$|\mathbf{w}|$"),
    ], fontsize=7, loc="lower left", ncol=2)
    fig.tight_layout()
    fig.savefig(path)


def tables(root, path):
    lines = ["% D2"]
    for method, label in (("blockCoupled", "monolithic"),
                          ("thetaNewmark", "segregated")):
        for co in ("0.5", "2"):
            case = root / "D2" / f"{method}_Co{co}"
            s = json.loads((case / "postProcessing/pmlVerification/"
                            "summary.json").read_text())
            st = s["setup"]
            centres = -st["domainHalfWidth"] + (np.arange(st["n"]) + 0.5) * st["h"]
            depth = np.clip(np.abs(centres) - st["halfWidth"], 0, None).max()
            sdt = st["sigmaMax"] * (depth / st["pmlThickness"]) ** 3 * st["deltaT"]
            if s["lastSampleTime"] < 0.049:
                final = (f"stopped, ${s['lastSampleTime'] * 1e3:.1f}"
                         "\\,\\mathrm{ms}$")
                cells = [final, "--", "--"]
            else:
                m = s["maxima"]
                cells = [sci(s["energy"]["final"]["Eac"]),
                         sci(m["psi"]["final"]), sci(m["w"]["final"])]
            lines.append(f"    {label} & ${co}$ & ${sdt:.2f}$ & "
                         f"{sci(s['residual']['maxAbsP'])} & "
                         + " & ".join(cells) + " \\\\")
    lines.append("% D3")
    for co in ("0.5", "2", "8"):
        for sigma in ("1.26e5", "2.53e5", "5.05e5"):
            m = json.loads((root / "D3" / f"monolithic_Co{co}_sigma{sigma}"
                            / "spectrum.json").read_text())
            j = json.loads((root / "D3" / f"monolithic_Co{co}_sigma{sigma}"
                            / "jordan.json").read_text())
            cells = [f"${co}$", f"${float(sigma) / 1e5:.2f}$",
                     f"${m['outerCellSigmaDt']:.3g}$",
                     sci(m["spectralRadius"] - 1, 1),
                     f"${j['nullity1']}/{j['nullity2']}$"]
            if sigma == "5.05e5":
                seg = json.loads((root / "D3" / f"segregated_Co{co}_sigma{sigma}"
                                  / "spectrum.json").read_text())
                imp = json.loads((root / "D3" / f"segregatedImplicit_Co{co}_sigma{sigma}"
                                  / "spectrum.json").read_text())
                cells += [f"${seg['spectralRadius']:.3f}$",
                          f"${seg['explicitDampingRadius']:.3f}$",
                          f"${imp['spectralRadius']:.3f}$"]
            else:
                cells += ["", "", ""]
            lines.append("    " + " & ".join(cells) + " \\\\")
    Path(path).write_text("\n".join(lines) + "\n")


def main():
    root = Path(sys.argv[1])
    figure(root, sys.argv[2])
    tables(root, sys.argv[3])


if __name__ == "__main__":
    main()
