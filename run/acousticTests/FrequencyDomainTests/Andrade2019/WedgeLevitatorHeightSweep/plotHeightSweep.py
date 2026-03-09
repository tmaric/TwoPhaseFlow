#!/usr/bin/env python3
import csv
import sys

import matplotlib.pyplot as plt


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: plotHeightSweep.py <input_csv> [output_png]", file=sys.stderr)
        return 1

    input_csv = sys.argv[1]
    output_png = sys.argv[2] if len(sys.argv) > 2 else "sweepResults/heightSweep_force.png"

    h_mm = []
    fn_ref_axisym = []
    with open(input_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            h_mm.append(float(row["height_mm"]))
            fn_ref_axisym.append(float(row["Fn_reflector1_axisym_N"]))

    if not h_mm:
        print(f"No data rows in {input_csv}", file=sys.stderr)
        return 1

    plt.figure(figsize=(7.0, 4.5))
    plt.plot(h_mm, fn_ref_axisym, "-o", markersize=3, linewidth=1)
    plt.xlabel("Height D [mm]")
    plt.ylabel("Reflector normal force Fn (axisym) [N]")
    plt.title("Wedge Levitator Height Sweep (frequency-domain)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_png, dpi=180)
    print(f"Wrote {output_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
