#!/usr/bin/env python3
"""Compare the resolved rigid-sphere force with the ideal Gorkov force."""

import argparse
import csv
from math import pi, sin
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--figure", required=True, type=Path)
    parser.add_argument("--frequency", required=True, type=float)
    parser.add_argument("--pressure-amplitude", required=True, type=float)
    parser.add_argument("--density", required=True, type=float)
    parser.add_argument("--sound-speed", required=True, type=float)
    parser.add_argument("--particle-y", required=True, type=float)
    parser.add_argument("--phase", default=0.0, type=float)
    return parser.parse_args()


def main():
    args = arguments()
    with args.input.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))

    radius = np.array([float(row["radius_m"]) for row in rows])
    numerical = np.array([float(row["numerical_force_y_N"]) for row in rows])
    k = 2.0*pi*args.frequency/args.sound_speed
    energy_density = (
        args.pressure_amplitude**2
        /(4.0*args.density*args.sound_speed**2)
    )

    # For an acoustically rigid sphere, f1=f2=1 and Phi=f1/3+f2/2=5/6.
    contrast_factor = 5.0/6.0
    analytical = (
        4.0*pi*radius**3*k*energy_density*contrast_factor
        * sin(2.0*(k*args.particle_y + args.phase))
    )
    relative_error = np.abs(numerical - analytical)/np.abs(analytical)

    fieldnames = list(rows[0]) + [
        "ka",
        "gorkov_force_y_N",
        "numerical_to_gorkov",
        "absolute_relative_error",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row, ka, exact, error, ratio in zip(
            rows, k*radius, analytical, relative_error, numerical/analytical
        ):
            row.update(
                {
                    "ka": f"{ka:.10g}",
                    "gorkov_force_y_N": f"{exact:.12g}",
                    "numerical_to_gorkov": f"{ratio:.10g}",
                    "absolute_relative_error": f"{error:.10g}",
                }
            )
            writer.writerow(row)

    fig, axes = plt.subplots(
        2, 1, figsize=(6.4, 6.2), sharex=True,
        gridspec_kw={"height_ratios": (2.0, 1.0)},
    )
    ka = k*radius
    axes[0].plot(ka, analytical*1e15, "k-", label="Gorkov prediction")
    axes[0].plot(ka, numerical*1e15, "o", color="#0072B2", label="Numerical")
    axes[0].set_ylabel(r"Axial force $F_y$ (fN)")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.25)

    axes[1].plot(ka, 100.0*relative_error, "o-", color="#D55E00")
    axes[1].set_xlabel(r"Size parameter $ka$")
    axes[1].set_ylabel("Relative error (%)")
    axes[1].grid(alpha=0.25)
    axes[1].set_ylim(bottom=0)

    fig.tight_layout()
    fig.savefig(args.figure, dpi=300)
    plt.close(fig)

    print(f"Wrote {args.output}")
    print(f"Wrote {args.figure}")
    print(f"Maximum relative error: {100.0*relative_error.max():.3f}%")


if __name__ == "__main__":
    main()
