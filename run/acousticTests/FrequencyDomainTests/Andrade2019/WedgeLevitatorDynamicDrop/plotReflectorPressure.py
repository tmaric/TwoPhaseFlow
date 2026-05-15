#!/usr/bin/env python3
from pathlib import Path
import sys

import matplotlib.pyplot as plt


def load_p0(path: Path):
    time = []
    p0 = []
    fn_axisym = []

    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split()
            time.append(float(cols[0]))
            fn_axisym.append(float(cols[2]))
            p0.append(float(cols[3]))

    return time, fn_axisym, p0


def main() -> int:
    case_dir = Path(__file__).resolve().parent
    input_file = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else case_dir / "postProcessing" / "reflectorPressure" / "p0.dat"
    )

    output_base = (
        Path(sys.argv[2])
        if len(sys.argv) > 2
        else case_dir / "postProcessing" / "reflectorPressure" / "p0"
    )

    if not input_file.is_file():
        print(f"Missing input file: {input_file}", file=sys.stderr)
        return 1

    time, fn_axisym, p0 = load_p0(input_file)

    if not time:
        print(f"No data rows in {input_file}", file=sys.stderr)
        return 1

    time_ms = [1.0e3*t for t in time]
    p0_kpa = [1.0e-3*p for p in p0]
    fn_mn = [1.0e3*f for f in fn_axisym]

    fig, axes = plt.subplots(2, 1, figsize=(6.5, 5.6), sharex=True)

    axes[0].plot(time_ms, p0_kpa, color="#1f4e79", linewidth=1.8)
    axes[0].set_ylabel("p0 [kPa]")
    axes[0].grid(True, alpha=0.25, linewidth=0.6)

    axes[1].plot(time_ms, fn_mn, color="#7a3b69", linewidth=1.8)
    axes[1].set_xlabel("time [ms]")
    axes[1].set_ylabel("F_R [mN]")
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
