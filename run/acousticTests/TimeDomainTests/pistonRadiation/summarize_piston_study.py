#!/usr/bin/env python3
"""Tabulate the piston study: on-axis and ring amplitudes against the
analytical baffled-piston solution (verify_piston_radiation.py output),
with the relative L2 error over the probe points and the mean ratio."""
import csv
import math
import sys
from pathlib import Path

import numpy as np


def errors(path):
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    sim, ana = data[:, 1], data[:, 2]
    return (math.sqrt(np.sum((sim - ana) ** 2) / np.sum(ana**2)),
            float(np.abs(sim - ana).max() / np.abs(ana).max()),
            float(np.median(sim / ana)), len(sim))


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "studies")
    rows = []
    for case in sorted(root.glob("monolithic_N*"), key=lambda c: int(c.name[12:])):
        if not (case / "verification_axis.csv").is_file():
            rows.append([case.name] + [""] * 8)
            continue
        a = errors(case / "verification_axis.csv")
        r = errors(case / "verification_ring.csv")
        rows.append([case.name[12:], f"{a[0]:.3e}", f"{a[1]:.3e}", f"{a[2]:.4f}",
                     a[3], f"{r[0]:.3e}", f"{r[1]:.3e}", f"{r[2]:.4f}", r[3]])
    header = ["cellsPerLambda", "axisRelL2", "axisRelMax", "axisMedianRatio",
              "axisPoints", "ringRelL2", "ringRelMax", "ringMedianRatio",
              "ringPoints"]
    with (root / "studyPiston.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    widths = [max(len(str(v)) for v in col) for col in zip(header, *rows)]
    for row in [header, *rows]:
        print("  ".join(str(v).rjust(w) for v, w in zip(row, widths)))


if __name__ == "__main__":
    main()
