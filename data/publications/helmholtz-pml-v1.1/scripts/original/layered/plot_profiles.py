#!/usr/bin/env python3
"""Pressure profiles of a solved layered case as one PDF figure.

Reads the comparison CSV written by postprocess_interface.py and the case's
studyProperties, and plots the requested panels (Pre, Pim, abs) of the
numerical and analytical pressure along the domain with the PML shaded.
"""
import argparse
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from postprocess_interface import parse_keyvals

COLUMNS = {"Pre": ("Pre_sim", "Pre_analytic", r"Re $P$ [Pa]"),
           "Pim": ("Pim_sim", "Pim_analytic", r"Im $P$ [Pa]"),
           "abs": ("p_abs_sim", "p_abs_analytic", r"$|P|$ [Pa]")}
HEIGHT = {1: 3.7, 2: 5.4, 3: 7.6}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", type=pathlib.Path)
    ap.add_argument("output", type=pathlib.Path)
    ap.add_argument("--panels", default="Pre,Pim", help="comma-separated subset of Pre, Pim, abs")
    args = ap.parse_args()
    panels = args.panels.split(",")
    props = parse_keyvals(args.case/"studyProperties")
    csv = sorted((args.case/"postProcessing"/"interfaceValidation").glob("*/pressureFieldComparison.csv"))
    if not csv:
        raise SystemExit(f"no pressureFieldComparison.csv in {args.case}")
    data = np.genfromtxt(csv[-1], delimiter=",", names=True)

    fig, axes = plt.subplots(len(panels), 1, sharex=True, figsize=(7.0, HEIGHT[len(panels)]), squeeze=False)
    for ax, panel in zip(axes[:, 0], panels):
        sim, ana, label = COLUMNS[panel]
        ax.axvspan(props["X_MAX"] - props["PML_L"], props["X_MAX"], color="0.9", lw=0)
        ax.plot(data["x_m"], data[ana], "k-", lw=2.0, label="Analytical")
        ax.plot(data["x_m"], data[sim], "--", color="C0", lw=1.8, label="Numerical")
        ax.set_xlim(props["X_MIN"], props["X_MAX"])
        ax.set_ylabel(label)
        ax.grid(alpha=0.25)
        ax.legend(loc="upper right", fontsize=8)
    axes[-1, 0].set_xlabel(r"$x$ [m]")
    fig.tight_layout()
    fig.savefig(args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
