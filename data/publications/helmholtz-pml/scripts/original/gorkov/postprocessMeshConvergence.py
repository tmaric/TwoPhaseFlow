#!/usr/bin/env python3
"""Postprocess the rigid-sphere near-body mesh refinement study."""

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
    parser.add_argument("--radius", required=True, type=float)
    parser.add_argument("--frequency", required=True, type=float)
    parser.add_argument("--pressure-amplitude", required=True, type=float)
    parser.add_argument("--density", required=True, type=float)
    parser.add_argument("--sound-speed", required=True, type=float)
    parser.add_argument("--particle-y", required=True, type=float)
    return parser.parse_args()


def main():
    args = arguments()
    with args.input.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))

    h_over_a = np.array([float(row["h_over_a"]) for row in rows])
    numerical = np.array([float(row["numerical_force_y_N"]) for row in rows])
    k = 2.0*pi*args.frequency/args.sound_speed
    energy_density = (
        args.pressure_amplitude**2
        /(4.0*args.density*args.sound_speed**2)
    )
    gorkov = (
        4.0*pi*args.radius**3*k*energy_density*(5.0/6.0)
        * sin(2.0*k*args.particle_y)
    )
    relative_difference = np.abs(numerical - gorkov)/abs(gorkov)
    change_from_finer = np.full(len(rows), np.nan)
    change_from_finer[:-1] = (
        np.abs(numerical[:-1] - numerical[1:])/np.abs(numerical[1:])
    )

    fieldnames = list(rows[0]) + [
        "gorkov_force_y_N",
        "difference_from_gorkov",
        "change_to_next_finer",
    ]
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for i, row in enumerate(rows):
            row.update(
                {
                    "gorkov_force_y_N": f"{gorkov:.12g}",
                    "difference_from_gorkov": f"{relative_difference[i]:.10g}",
                    "change_to_next_finer": (
                        "" if np.isnan(change_from_finer[i])
                        else f"{change_from_finer[i]:.10g}"
                    ),
                }
            )
            writer.writerow(row)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))
    force_scale = 1e15
    axes[0].axhline(gorkov*force_scale, color="black", label="Gorkov")
    axes[0].plot(
        h_over_a, numerical*force_scale, "o-", color="#0072B2",
        label="Numerical",
    )
    axes[0].invert_xaxis()
    axes[0].set_xlabel(r"Near-sphere resolution $h_s/a$")
    axes[0].set_ylabel(r"Axial force $F_y$ (fN)")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)

    axes[1].semilogy(
        h_over_a, 100.0*relative_difference, "o-", color="#D55E00"
    )
    axes[1].invert_xaxis()
    axes[1].set_xlabel(r"Near-sphere resolution $h_s/a$")
    axes[1].set_ylabel("Difference from Gorkov (%)")
    axes[1].grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(args.figure, dpi=300)
    plt.close(fig)

    print(f"Wrote {args.output}")
    print(f"Wrote {args.figure}")


if __name__ == "__main__":
    main()
