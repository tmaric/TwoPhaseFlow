#!/usr/bin/env python3
from pathlib import Path
import sys

import matplotlib.pyplot as plt


def load_aspect_ratio(path: Path):
    time = []
    h_over_v = []
    v_over_h = []
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
            h_over_v.append(float(cols[3]))
            v_over_h.append(float(cols[4]))

    return time, h_over_v, v_over_h, horizontal_diameter, vertical_diameter


def load_volume(path: Path):
    time = []
    volume = []
    equivalent_radius = []

    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split()
            time.append(float(cols[0]))
            volume.append(float(cols[1]))
            equivalent_radius.append(float(cols[3]))

    return time, volume, equivalent_radius


def main() -> int:
    case_dir = Path(__file__).resolve().parent
    input_file = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else case_dir / "postProcessing" / "dropAspectRatio" / "aspectRatio.dat"
    )

    output_base = (
        Path(sys.argv[2])
        if len(sys.argv) > 2
        else case_dir / "postProcessing" / "dropAspectRatio" / "aspectRatio"
    )
    volume_file = (
        Path(sys.argv[3])
        if len(sys.argv) > 3
        else case_dir / "postProcessing" / "dropVolume" / "volume.dat"
    )

    if not input_file.is_file():
        print(f"Missing input file: {input_file}", file=sys.stderr)
        return 1
    if not volume_file.is_file():
        print(f"Missing volume file: {volume_file}", file=sys.stderr)
        return 1

    time, h_over_v, v_over_h, horizontal_diameter, vertical_diameter = (
        load_aspect_ratio(input_file)
    )
    volume_time, volume, equivalent_radius = load_volume(volume_file)

    if not time:
        print(f"No data rows in {input_file}", file=sys.stderr)
        return 1
    if not volume_time:
        print(f"No data rows in {volume_file}", file=sys.stderr)
        return 1

    horizontal_mm = [1.0e3*d for d in horizontal_diameter]
    vertical_mm = [1.0e3*d for d in vertical_diameter]
    volume_mm3 = [1.0e9*v for v in volume]
    equivalent_radius_mm = [1.0e3*r for r in equivalent_radius]

    fig, axes = plt.subplots(3, 1, figsize=(6.5, 7.2), sharex=True)

    axes[0].plot(time, h_over_v, color="#1f4e79", linewidth=1.8)
    axes[0].plot(time, v_over_h, color="#b23a48", linewidth=1.4, linestyle="--")
    axes[0].set_ylabel("aspect ratio [-]")
    axes[0].legend(["H/V", "V/H"], frameon=False)
    axes[0].grid(True, alpha=0.25, linewidth=0.6)

    axes[1].plot(time, horizontal_mm, color="#1f4e79", linewidth=1.8)
    axes[1].plot(time, vertical_mm, color="#b23a48", linewidth=1.4, linestyle="--")
    axes[1].set_ylabel("diameter [mm]")
    axes[1].legend(["horizontal", "vertical"], frameon=False)
    axes[1].grid(True, alpha=0.25, linewidth=0.6)

    axes[2].plot(volume_time, volume_mm3, color="#2f6f4e", linewidth=1.8)
    axes_radius = axes[2].twinx()
    axes_radius.plot(
        volume_time,
        equivalent_radius_mm,
        color="#7a3b69",
        linewidth=1.4,
        linestyle="--",
    )
    axes[2].set_xlabel("pseudo-time [-]")
    axes[2].set_ylabel("volume [mm3]")
    axes_radius.set_ylabel("equiv. radius [mm]")
    axes[2].grid(True, alpha=0.25, linewidth=0.6)

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
