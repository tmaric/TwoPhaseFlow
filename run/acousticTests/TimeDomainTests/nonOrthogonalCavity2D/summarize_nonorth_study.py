#!/usr/bin/env python3
"""Tabulate the non-orthogonal cavity studies written by runStudy."""
import csv
import json
import math
import re
import sys
from pathlib import Path


def summary(case):
    s = json.loads((case / "postProcessing/cavity2dVerification/summary.json").read_text())
    text = (case / "log.checkMesh").read_text()
    m = re.search(r"Mesh non-orthogonality Max: ([0-9.eE+-]+) average: ([0-9.eE+-]+)", text)
    s["nonOrth"] = (float(m.group(1)), float(m.group(2))) if m else (0.0, 0.0)
    return s


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
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "studies")
    if (root / "A").is_dir():
        rows = []
        for eps in sorted({c.name.split("_")[0] for c in (root / "A").iterdir()},
                          key=lambda e: float(e[3:])):
            items = sorted((summary(c) for c in (root / "A").glob(f"{eps}_N*")
                            if (c / "postProcessing").is_dir()),
                           key=lambda s: s["setup"]["cells"])
            prev = None
            for s in items:
                e = s["E2"]
                rows.append([
                    eps[3:], s["setup"]["cells"], f"{s['nonOrth'][0]:.1f}",
                    f"{s['nonOrth'][1]:.1f}", f"{e['pressure']:.3e}",
                    f"{math.log(prev['pressure'] / e['pressure'], 2):.2f}" if prev else "",
                    f"{e['velocity']:.3e}",
                    f"{math.log(prev['velocity'] / e['velocity'], 2):.2f}" if prev else "",
                    f"{s['modal']['relativeFrequencyError']:+.3e}",
                    f"{s['energyDrift']:+.2e}"])
                prev = e
        write(root / "studyA.csv",
              ["distortion", "cells", "nonOrthMax", "nonOrthAvg", "E2p", "q",
               "E2U", "q", "freqError", "energyDrift"], rows)
    if (root / "B").is_dir():
        rows = []
        for c in sorted((root / "B").iterdir()):
            if not (c / "postProcessing/cavity2dVerification").is_dir():
                rows.append([c.name, "", "", "", "no summary", "", ""])
                continue
            s = summary(c)
            rows.append([c.name, f"{s['setup']['courant']:.2f}",
                         f"{s['lastTime']:.4g}", s["finite"],
                         f"{s['modal']['growthPerPeriod'] - 1:+.3e}",
                         f"{s['energyDrift']:+.3e}",
                         f"{s['E2']['pressureLastPeriod']:.3e}"])
        write(root / "studyB.csv",
              ["case", "Co", "tLast", "finite", "modalGrowthPerPeriod",
               "energyDrift", "E2pLastPeriod"], rows)


if __name__ == "__main__":
    main()
