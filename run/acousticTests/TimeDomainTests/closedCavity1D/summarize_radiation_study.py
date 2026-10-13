#!/usr/bin/env python3
"""Tabulate the cavity radiation-average studies written by runRadiationStudy."""
import csv
import json
import math
import sys
from pathlib import Path


def summary(case):
    return json.loads((case / "postProcessing/radiationVerification/"
                       "summary.json").read_text())


def fmt(x):
    return f"{x:.3e}"


def order(coarse, fine):
    return f"{math.log(coarse / fine) / math.log(2.0):.2f}"


def write(path, header, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"\n{path}")
    widths = [max(len(str(v)) for v in col) for col in zip(header, *rows)]
    for row in [header, *rows]:
        print("  ".join(str(v).rjust(w) for v, w in zip(row, widths)))


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "studies/radiation")
    if (root / "A").is_dir():
        items = sorted((summary(c) for c in (root / "A").iterdir()),
                       key=lambda s: s["setup"]["radiationPeriods"])
        rows = [[s["setup"]["radiationPeriods"], s["windowSamples"],
                 f"{s['windowDiscretePeriods']:.4f}",
                 fmt(s["biasBoundDiscrete"]),
                 fmt(s["errors"]["discrete"]["pr"]),
                 fmt(s["errors"]["exact"]["pr"]),
                 fmt(s["errors"]["mean"]["pr"]),
                 fmt(s["errors"]["mean"]["Rxx"])] for s in items]
        write(root / "studyA.csv",
              ["M0", "N_s", "M_discrete", "biasBound", "pr_vs_discrete",
               "pr_vs_exact", "pr_vs_mean", "Rxx_vs_mean"], rows)
    if (root / "B").is_dir():
        items = sorted((summary(c) for c in (root / "B").iterdir()),
                       key=lambda s: s["setup"]["cellsPerWavelength"])
        rows = []
        for i, s in enumerate(items):
            e, prev = s["errors"], items[i - 1]["errors"] if i else None
            rows.append([
                s["setup"]["cellsPerWavelength"], s["setup"]["stepsPerPeriod"],
                fmt(s["biasBoundDiscrete"]), fmt(e["discrete"]["pr"]),
                fmt(e["exact"]["pr"]),
                order(prev["exact"]["pr"], e["exact"]["pr"]) if i else "",
                fmt(e["mean"]["pr"]),
                order(prev["mean"]["pr"], e["mean"]["pr"]) if i else "",
                fmt(e["mean"]["Rxx"]),
                order(prev["mean"]["Rxx"], e["mean"]["Rxx"]) if i else "",
            ])
        write(root / "studyB.csv",
              ["N_lambda", "N_T", "biasBound", "pr_vs_discrete", "pr_vs_exact",
               "order", "pr_vs_mean", "order", "Rxx_vs_mean", "order"], rows)


if __name__ == "__main__":
    main()
