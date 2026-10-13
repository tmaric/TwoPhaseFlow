#!/usr/bin/env python3
"""Tabulate the two-layer radiation-average study written by runRadiationStudy."""
import csv
import json
import math
import sys
from pathlib import Path

GROUPS = ("gasInterior", "gasInterfaceCell", "liquidInterfaceCell",
          "liquidInterior")


def summary(case):
    return json.loads((case / "postProcessing/radiationVerification/"
                       "summary.json").read_text())


def order(coarse, fine):
    return f"{math.log(coarse / fine) / math.log(2.0):.2f}"


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "studies/radiation")
    header = ["mode", "cellsPerUnit", "N_T", "discreteMax"]
    for g in GROUPS:
        header += [f"pr_{g}", "order", f"Rxx_{g}", "order"]
    rows = []
    for mode in ("lower", "upper"):
        items = sorted((summary(c) for c in (root / "C").glob(f"{mode}_*")),
                       key=lambda s: s["setup"]["cellsPerUnit"])
        for i, s in enumerate(items):
            e, prev = s["errors"]["mean"], items[i - 1]["errors"]["mean"] if i else None
            disc = max(v for g in s["errors"]["discrete"].values()
                       for v in g.values())
            row = [mode, s["setup"]["cellsPerUnit"],
                   s["setup"]["stepsPerPeriod"], f"{disc:.1e}"]
            for g in GROUPS:
                for q in ("pr", "Rxx"):
                    row += [f"{e[g][q]:.3e}",
                            order(prev[g][q], e[g][q]) if i else ""]
            rows.append(row)
    path = root / "studyC.csv"
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    print(path)
    widths = [max(len(str(v)) for v in col) for col in zip(header, *rows)]
    for row in [header, *rows]:
        print("  ".join(str(v).rjust(w) for v, w in zip(row, widths)))


if __name__ == "__main__":
    main()
