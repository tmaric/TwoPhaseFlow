#!/usr/bin/env python3
"""Compare empty-field Gorkov force with hard-hole integrated force."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"{path} has no header")
        return list(reader)


def as_float(row: dict[str, str], key: str, path: Path) -> float:
    try:
        return float(row[key])
    except KeyError as exc:
        raise SystemExit(f"{path} has no column {key!r}") from exc


def main() -> int:
    empty_summary = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else Path("emptyGorkovStudy/summary.tsv")
    )
    hole_summary = (
        Path(sys.argv[2])
        if len(sys.argv) > 2
        else Path("../WedgeLevitatorGorkovHole/solidBallGorkovStudy/summary.tsv")
    )
    output_png = (
        Path(sys.argv[3])
        if len(sys.argv) > 3
        else empty_summary.with_name("gorkovVsHoleForce.png")
    )
    output_tsv = (
        Path(sys.argv[4])
        if len(sys.argv) > 4
        else empty_summary.with_name("gorkovVsHoleForce.tsv")
    )

    empty_rows = read_table(empty_summary)
    hole_rows = read_table(hole_summary)

    empty_by_radius = {
        as_float(row, "radius_um", empty_summary): as_float(
            row, "empty_Fy_gorkov", empty_summary
        )
        for row in empty_rows
    }
    hole_by_radius = {
        as_float(row, "radius_um", hole_summary): as_float(
            row, "sim_Fy_axisym", hole_summary
        )
        for row in hole_rows
    }

    radii = sorted(set(empty_by_radius) & set(hole_by_radius))
    if not radii:
        raise SystemExit("No common radius values between the two summaries")

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["radius_um", "empty_Fy_gorkov", "hole_Fy_axisym"])
        for radius in radii:
            writer.writerow([radius, empty_by_radius[radius], hole_by_radius[radius]])

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.6, 4.3), constrained_layout=True)
    ax.plot(
        radii,
        [empty_by_radius[radius] for radius in radii],
        marker="o",
        linewidth=1.8,
        label="Empty: Gorkov point force",
    )
    ax.plot(
        radii,
        [hole_by_radius[radius] for radius in radii],
        marker="s",
        linewidth=1.8,
        label="GorkovHole: integrated pressure force",
    )
    ax.axhline(0.0, color="0.35", linewidth=0.8)
    ax.set_xlabel("Sphere radius [um]")
    ax.set_ylabel("Acoustic radiation force Fy [N]")
    ax.set_title("Gorkov point force vs hard-hole integrated force")
    ax.grid(True, color="0.85", linewidth=0.8)
    ax.legend()

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=200)
    print(f"Wrote {output_png}")
    print(f"Wrote {output_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
