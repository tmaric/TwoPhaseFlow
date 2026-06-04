#!/usr/bin/env python3
from pathlib import Path
import sys

import matplotlib.pyplot as plt


def load_shape(path: Path):
    time = []
    horizontal_diameter = []
    vertical_diameter = []

    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split()
            time.append(float(cols[0]))
            horizontal_diameter.append(float(cols[1]))
            vertical_diameter.append(float(cols[2]))

    return time, horizontal_diameter, vertical_diameter


def load_position(path: Path):
    time = []
    center_y = []

    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split()
            time.append(float(cols[0]))
            center_y.append(float(cols[1]))

    return time, center_y


def main() -> int:
    case_dir = Path(__file__).resolve().parent
    shape_file = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else case_dir / "postProcessing" / "dropAspectRatio" / "aspectRatio.dat"
    )
    position_file = (
        Path(sys.argv[2])
        if len(sys.argv) > 2
        else case_dir / "postProcessing" / "dropCenter" / "center.dat"
    )
    output_png = (
        Path(sys.argv[3])
        if len(sys.argv) > 3
        else case_dir / "postProcessing" / "shapePosition.png"
    )

    if not shape_file.is_file():
        print(f"Missing input file: {shape_file}", file=sys.stderr)
        return 1
    if not position_file.is_file():
        print(f"Missing input file: {position_file}", file=sys.stderr)
        return 1

    shape_time, horizontal_diameter, vertical_diameter = load_shape(shape_file)
    position_time, center_y = load_position(position_file)

    if not shape_time:
        print(f"No data rows in {shape_file}", file=sys.stderr)
        return 1
    if not position_time:
        print(f"No data rows in {position_file}", file=sys.stderr)
        return 1

    horizontal_mm = [1.0e3*d for d in horizontal_diameter]
    vertical_mm = [1.0e3*d for d in vertical_diameter]
    center_y_mm = [1.0e3*y for y in center_y]

    fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.4), sharex=False)

    axes[0].plot(shape_time, horizontal_mm, color="#1f4e79", linewidth=1.8)
    axes[0].plot(
        shape_time,
        vertical_mm,
        color="#b23a48",
        linewidth=1.5,
        linestyle="--",
    )
    axes[0].set_ylabel("diameter [mm]")
    axes[0].legend(["horizontal", "vertical"], frameon=False)
    axes[0].grid(True, alpha=0.25, linewidth=0.6)

    axes[1].plot(position_time, center_y_mm, color="#2f6f4e", linewidth=1.8)
    axes[1].set_xlabel("pseudo-time [-]")
    axes[1].set_ylabel("centroid y [mm]")
    axes[1].grid(True, alpha=0.25, linewidth=0.6)

    fig.tight_layout()

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=220)

    print(f"Wrote {output_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
