#!/usr/bin/env python3
"""Compare resolved and Gorkov forces as a rigid sphere moves vertically."""

import argparse
import csv
from math import pi
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--figure", required=True, type=Path)
    parser.add_argument("--radius", required=True, type=float)
    parser.add_argument("--frequency", required=True, type=float)
    parser.add_argument("--pressure-amplitude", required=True, type=float)
    parser.add_argument("--density", required=True, type=float)
    parser.add_argument("--sound-speed", required=True, type=float)
    parser.add_argument("--height", required=True, type=float)
    return parser.parse_args()


def main():
    args = arguments()
    with args.input.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))

    y = np.array([float(row["position_m"]) for row in rows])
    numerical = np.array([float(row["numerical_force_y_N"]) for row in rows])
    k = 2.0*pi*args.frequency/args.sound_speed
    energy_density = (
        args.pressure_amplitude**2
        /(4.0*args.density*args.sound_speed**2)
    )
    peak_force = 4.0*pi*args.radius**3*k*energy_density*(5.0/6.0)
    analytical = peak_force*np.sin(2.0*k*y)
    normalized_difference = np.abs(numerical - analytical)/abs(peak_force)

    fieldnames = list(rows[0]) + [
        "gorkov_force_y_N",
        "difference_normalized_by_peak_force",
    ]
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row, exact, difference in zip(
            rows, analytical, normalized_difference
        ):
            row.update(
                {
                    "gorkov_force_y_N": f"{exact:.12g}",
                    "difference_normalized_by_peak_force": f"{difference:.10g}",
                }
            )
            writer.writerow(row)

    y_line = np.linspace(0.0, args.height, 500)
    fig, axes = plt.subplots(
        2, 1, figsize=(6.4, 6.2), sharex=True,
        gridspec_kw={"height_ratios": (2.0, 1.0)},
    )
    axes[0].plot(
        1e3*y_line, 1e15*peak_force*np.sin(2.0*k*y_line),
        "k-", label="Gorkov prediction",
    )
    axes[0].plot(
        1e3*y, 1e15*numerical, "o", color="#0072B2", label="Numerical"
    )
    axes[0].axhline(0, color="0.5", linewidth=0.8)
    axes[0].set_ylabel(r"Axial force $F_y$ (fN)")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.25)

    axes[1].plot(
        1e3*y, 100.0*normalized_difference, "o-", color="#D55E00"
    )
    axes[1].set_xlabel(r"Sphere-centre position $y$ (mm)")
    axes[1].set_ylabel(r"$|F_h-F_G|/|F_G|_{\max}$ (%)")
    axes[1].set_ylim(bottom=0)
    axes[1].grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(args.figure, dpi=300)
    plt.close(fig)

    rms = np.sqrt(np.mean(normalized_difference**2))
    print(f"Wrote {args.output}")
    print(f"Wrote {args.figure}")
    print(f"Maximum peak-normalized difference: {100*normalized_difference.max():.3f}%")
    print(f"RMS peak-normalized difference: {100*rms:.3f}%")


if __name__ == "__main__":
    main()
