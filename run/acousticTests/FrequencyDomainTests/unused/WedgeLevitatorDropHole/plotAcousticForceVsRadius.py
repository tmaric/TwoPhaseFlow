#!/usr/bin/env python3
"""Plot acoustic radiation force versus sphere radius for GorkovHole sweeps."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


def read_summary(path: Path) -> tuple[list[float], list[float]]:
    radii_um: list[float] = []
    force_y: list[float] = []

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"radius_um", "sim_Fy_axisym"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            fields = ", ".join(reader.fieldnames or [])
            raise SystemExit(
                f"{path} must contain tab-separated columns "
                f"{', '.join(sorted(required))}; found: {fields}"
            )

        for row in reader:
            radii_um.append(float(row["radius_um"]))
            force_y.append(float(row["sim_Fy_axisym"]))

    if not radii_um:
        raise SystemExit(f"{path} has no data rows")

    return radii_um, force_y


def main() -> int:
    summary = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("solidBallGorkovStudy/summary.tsv")
    output = (
        Path(sys.argv[2])
        if len(sys.argv) > 2
        else summary.with_name("acousticForceVsRadius.png")
    )

    radii_um, force_y = read_summary(summary)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    ax.plot(radii_um, force_y, marker="o", linewidth=1.8)
    ax.axhline(0.0, color="0.35", linewidth=0.8)
    ax.set_xlabel("Sphere radius [um]")
    ax.set_ylabel("Axisymmetric acoustic force Fy [N]")
    ax.set_title("Acoustic radiation force vs sphere radius")
    ax.grid(True, color="0.85", linewidth=0.8)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200)
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
