#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    source = pathlib.Path(sys.argv[1])
    out = pathlib.Path(sys.argv[2])
    with source.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    groups: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(row["offset"], []).append(row)
    for group in groups.values():
        group.sort(key=lambda row: float(row["NX"]))
        for index, row in enumerate(group):
            if index == 0:
                row["pressureOrder"] = ""
            else:
                previous = group[index-1]
                row["pressureOrder"] = f"{math.log(float(previous['pressureRelL2'])/float(row['pressureRelL2']))/math.log(float(row['NX'])/float(previous['NX'])):.8g}"
    fieldnames = [key for key in rows[0] if key != "pressureOrder"] + ["pressureOrder"]
    with source.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(groups, key=float):
            writer.writerows(groups[key])

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    for offset, group in sorted(groups.items(), key=lambda item: float(item[0])):
        h = [float(row["dx"]) for row in group]
        axes[0].loglog(h, [float(row["pressureRelL2"]) for row in group], "o-", label=rf"$\delta={offset}$")
        axes[1].loglog(h, [float(row["RAbsError"]) for row in group], "o-", label=rf"$\delta={offset}$")
    for ax, ylabel in zip(axes, [r"$E_2(P)$", r"$|R_h-R|$"]):
        ax.set_xlabel(r"$\Delta x$ [m]")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.3)
        ax.invert_xaxis()
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out / "interface_convergence.png", dpi=220)

    with (out / "interface_convergence_table.tex").open("w", encoding="utf-8") as stream:
        stream.write("\\begin{table}[htbp]\n\\centering\n\\begin{tabular}{rrrrrr}\n\\toprule\n")
        stream.write("$\delta$ & $N_x$ & $E_2(P)$ & $q_P$ & $|R_h-R|$ & $|T_h-T|$ \\\\\n\\midrule\n")
        for offset, group in sorted(groups.items(), key=lambda item: float(item[0])):
            for row in group[-4:]:
                q = "--" if not row["pressureOrder"] else f"{float(row['pressureOrder']):.2f}"
                stream.write(f"{float(offset):.2f} & {float(row['NX']):.0f} & {float(row['pressureRelL2']):.3e} & {q} & {float(row['RAbsError']):.3e} & {float(row['TAbsError']):.3e} \\\\\n")
        stream.write("\\bottomrule\n\\end{tabular}\n")
        stream.write("\\caption{Layered-interface convergence for interface offsets $\delta=0$, $0.25$, and $0.5$ cell widths from a mesh face.}\n")
        stream.write("\\label{tab:interfaceConvergence}\n\\end{table}\n")


if __name__ == "__main__":
    main()
