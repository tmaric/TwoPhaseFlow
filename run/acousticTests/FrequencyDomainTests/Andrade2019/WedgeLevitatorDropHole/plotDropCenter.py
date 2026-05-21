#!/usr/bin/env python3
from pathlib import Path
import sys

import matplotlib.pyplot as plt


def load_center(path: Path):
    time = []
    center_y = []
    vertical_diameter = []

    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split()
            time.append(float(cols[0]))
            center_y.append(float(cols[1]))
            vertical_diameter.append(float(cols[5]))

    return time, center_y, vertical_diameter


def main() -> int:
    case_dir = Path(__file__).resolve().parent
    input_file = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else case_dir / "postProcessing" / "dropCenter" / "center.dat"
    )

    output_base = (
        Path(sys.argv[2])
        if len(sys.argv) > 2
        else case_dir / "postProcessing" / "dropCenter" / "center"
    )

    if not input_file.is_file():
        print(f"Missing input file: {input_file}", file=sys.stderr)
        return 1

    time, center_y, vertical_diameter = load_center(input_file)

    if not time:
        print(f"No data rows in {input_file}", file=sys.stderr)
        return 1

    center_y_mm = [1.0e3*y for y in center_y]
    vertical_diameter_mm = [1.0e3*d for d in vertical_diameter]

    fig, axes = plt.subplots(2, 1, figsize=(6.5, 5.4), sharex=True)

    axes[0].plot(time, center_y_mm, color="#1f4e79", linewidth=1.8)
    axes[0].set_ylabel("centerY [mm]")
    axes[0].grid(True, alpha=0.25, linewidth=0.6)

    axes[1].plot(time, vertical_diameter_mm, color="#b23a48", linewidth=1.6)
    axes[1].set_xlabel("pseudo-time [-]")
    axes[1].set_ylabel("vertical diameter [mm]")
    axes[1].grid(True, alpha=0.25, linewidth=0.6)

    fig.tight_layout()

    output_base.parent.mkdir(parents=True, exist_ok=True)
    output_png = output_base.with_suffix(".png")
    output_pdf = output_base.with_suffix(".pdf")
    fig.savefig(output_png, dpi=220)
    fig.savefig(output_pdf)

    print(f"Wrote {output_png}")
    print(f"Wrote {output_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
