#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


def order(previous: dict[str, str], current: dict[str, str], key: str) -> float:
    e0, e1 = float(previous[key]), float(current[key])
    n0 = float(previous["cellsPerWavelength"])
    n1 = float(current["cellsPerWavelength"])
    return math.log(e0/e1)/math.log(n1/n0)


def main() -> None:
    csv_path = pathlib.Path(sys.argv[1])
    out_dir = pathlib.Path(sys.argv[2])
    with csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    groups: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault((row["meshFamily"], row["boundaryMode"]), []).append(row)
    for group in groups.values():
        group.sort(key=lambda row: float(row["cellsPerWavelength"]))
        for index, row in enumerate(group):
            row["pressureOrder"] = "" if index == 0 else f"{order(group[index-1], row, 'pressureRelL2'):.8g}"
            row["velocityOrder"] = "" if index == 0 else f"{order(group[index-1], row, 'velocityRelL2'):.8g}"

    order_fields = ["pressureOrder", "velocityOrder"]
    fieldnames = [key for key in rows[0] if key not in order_fields] + order_fields
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(groups):
            writer.writerows(groups[key])

    colors = {
        "orthogonal": "#0072B2",
        "warpedInterior": "#009E73",
    }
    plot_family_labels = {
        "orthogonal": "Orthogonal",
        "warpedInterior": "Non-orthogonal",
    }
    table_family_labels = {
        "orthogonal": "Orthogonal",
        "warpedInterior": "Non-orthogonal interior",
    }
    boundary_labels = {"dirichlet": "Dirichlet", "mixed": "Mixed"}
    markers = {"dirichlet": "o", "mixed": "s"}
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    for (family, boundary), group in sorted(groups.items()):
        x = [1/float(row["cellsPerWavelength"]) for row in group]
        label = f"{plot_family_labels[family]}, {boundary_labels[boundary]}"
        axes[0].loglog(x, [float(row["pressureRelL2"]) for row in group], marker=markers[boundary], color=colors[family], ls="-" if boundary == "dirichlet" else "--", label=label)
        axes[1].loglog(x, [float(row["velocityRelL2"]) for row in group], marker=markers[boundary], color=colors[family], ls="-" if boundary == "dirichlet" else "--", label=label)
    for ax, ylabel in zip(axes, [r"$E_2(P)$", r"$E_2(\mathbf{u})$"]):
        ax.set_xlabel(r"$h/\lambda$")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.3)
        ax.invert_xaxis()
        tick_cells = (16, 32, 64, 96)
        ax.set_xticks([1/value for value in tick_cells])
        ax.set_xticklabels([rf"$1/{value}$" for value in tick_cells])
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    axes[0].legend(
        loc="upper right",
        fontsize=7,
        frameon=True,
        framealpha=1.0,
        edgecolor="0.75",
        borderpad=0.25,
        labelspacing=0.25,
        handlelength=1.6,
        handletextpad=0.4,
    )
    fig.tight_layout()
    fig.savefig(out_dir / "homogeneous_convergence.png", dpi=220)

    selected = []
    for key in sorted(groups):
        group = groups[key]
        selected.extend(group[-4:])
    with (out_dir / "homogeneous_convergence_table.tex").open("w", encoding="utf-8") as stream:
        stream.write("\\begin{table}[htbp]\n\\centering\n")
        stream.write("\\small\n\\setlength{\\tabcolsep}{4pt}\n")
        stream.write("\\begin{tabular}{llrrrrr}\n\\toprule\n")
        stream.write("Mesh & BC & $N_\\lambda$ & $E_2(P)$ & $q_P$ & $E_2(\\mathbf{u})$ & $q_{\\mathbf{u}}$ \\\\\n\\midrule\n")
        for row in selected:
            stream.write(
                f"{table_family_labels[row['meshFamily']]} & "
                f"{boundary_labels[row['boundaryMode']]} & "
                f"{float(row['cellsPerWavelength']):.0f} & "
                f"{float(row['pressureRelL2']):.3e} & "
                f"{float(row['pressureOrder']):.2f} & "
                f"{float(row['velocityRelL2']):.3e} & "
                f"{float(row['velocityOrder']):.2f} \\\\\n"
            )
        stream.write("\\bottomrule\n\\end{tabular}\n")
        stream.write("\\caption{Homogeneous plane-wave convergence. $N_\\lambda$ denotes cells per wavelength; $q_P$ and $q_{\\mathbf{u}}$ are the observed pressure and velocity orders.}\n")
        stream.write("\\label{tab:homogeneousConvergence}\n\\end{table}\n")


if __name__ == "__main__":
    main()
