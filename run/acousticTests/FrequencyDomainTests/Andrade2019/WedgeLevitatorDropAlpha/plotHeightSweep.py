#!/usr/bin/env python3
import csv
from pathlib import Path
import sys

import matplotlib.pyplot as plt


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: plotHeightSweep.py <input_csv> [output_base]", file=sys.stderr)
        return 1

    input_csv = sys.argv[1]
    output_base = Path(sys.argv[2] if len(sys.argv) > 2 else "sweepResults/reproducedFig4")

    series = {}
    with open(input_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = row["series_key"]
            bucket = series.setdefault(
                key,
                {
                    "label": row["series_label"],
                    "h_mm": [],
                    "fn_mN": [],
                },
            )
            bucket["h_mm"].append(float(row["height_mm"]))
            bucket["fn_mN"].append(1.0e3 * float(row["Fn_reflector1_axisym_N"]))

    if not series:
        print(f"No data rows in {input_csv}", file=sys.stderr)
        return 1

    fig, ax = plt.subplots(figsize=(6.0, 4.0))

    style_map = {
        "empty": {"color": "black", "linestyle": "-", "linewidth": 2.0, "label": "No object"},
        "ab1": {"color": "#2f6df6", "linestyle": "-", "linewidth": 1.6, "label": "Sphere (a/b=1)"},
        "ab2": {"color": "#44a35f", "linestyle": "-", "linewidth": 1.6, "label": "Oblate spheroid (a/b=2)"},
        "ab3": {"color": "#ff4aa2", "linestyle": "-", "linewidth": 1.6, "label": "Oblate spheroid (a/b=3)"},
    }

    for key in sorted(series):
        bucket = series[key]
        paired = sorted(zip(bucket["h_mm"], bucket["fn_mN"]))
        h_mm = [p[0] for p in paired]
        fn_mN = [p[1] for p in paired]
        style = dict(style_map.get(key, {"linewidth": 1.4, "label": bucket["label"]}))
        label = style.pop("label")
        ax.plot(h_mm, fn_mN, label=label, **style)

    ax.set_xlabel("H [mm]")
    ax.set_ylabel("Radiation force on reflector [mN]")
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.legend(frameon=False, fontsize=9)
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
