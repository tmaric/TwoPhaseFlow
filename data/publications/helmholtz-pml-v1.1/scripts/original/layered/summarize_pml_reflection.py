#!/usr/bin/env python3
from __future__ import annotations

import csv
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
    groups = {name: [row for row in rows if row["study"] == name] for name in ("damping", "thickness", "resolution", "order")}

    fig, axes = plt.subplots(2, 2, figsize=(8.8, 6.8))
    labels = {
        "damping": r"$\sigma_{\max}/\omega$",
        "thickness": r"$L_{\mathrm{PML}}/\lambda$",
        "resolution": r"cells per wavelength",
        "order": r"polynomial order",
    }
    for ax, name in zip(axes.flat, groups):
        group = groups[name]
        x = [float(row["value"]) for row in group]
        y = [float(row["reflectionMagnitude"]) for row in group]
        ax.semilogy(x, y, "o-", color="#0072B2")
        ax.set_xlabel(labels[name])
        ax.set_ylabel(r"$|B/A|$")
        ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "pml_reflection_sensitivity.png", dpi=220)

    with (out / "pml_reflection_table.tex").open("w", encoding="utf-8") as stream:
        stream.write("\\begin{table}[htbp]\n\\centering\n\\begin{tabular}{lrrr}\n\\toprule\n")
        stream.write("Study & Value & $|B/A|$ & fit residual \\\\\n\\midrule\n")
        for name in groups:
            for row in groups[name]:
                stream.write(f"{name.capitalize()} & {float(row['value']):.3g} & {float(row['reflectionMagnitude']):.3e} & {float(row['fitResidual']):.3e} \\\\\n")
        stream.write("\\bottomrule\n\\end{tabular}\n")
        stream.write("\\caption{Normal-incidence PML reflection obtained by fitting forward and backward waves in the physical domain.}\n")
        stream.write("\\label{tab:pmlReflection}\n\\end{table}\n")


if __name__ == "__main__":
    main()
