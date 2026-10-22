#!/usr/bin/env python3
"""Tabulate the wall-layer study written by runWallLayerStudy."""
import csv
import json
import math
import re
import sys
from pathlib import Path


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "studies/wallLayer")
    rows = []
    for scheme in ("orthogonal", "monolithic", "segregated"):
        for eps in ("0.15", "0.3", "0.6"):
            prev = None
            for n in (20, 40, 80, 160):
                case = root / f"{scheme}_eps{eps}_N{n}"
                path = case / "postProcessing/cavity2dVerification/summary.json"
                if not path.is_file():
                    rows.append([scheme, eps, n, "", "", "", "", "", "no summary"])
                    continue
                s = json.loads(path.read_text())
                m = re.search(r"Mesh non-orthogonality Max: ([0-9.eE+-]+) average: ([0-9.eE+-]+)",
                              (case / "log.checkMesh").read_text())
                e = s["E2"]["pressure"]
                ok = math.isfinite(e) and e < 10
                rows.append([scheme, eps, n, f"{float(m.group(1)):.1f}",
                             f"{float(m.group(2)):.1f}", f"{e:.3e}",
                             f"{math.log(prev / e, 2):.2f}" if prev and ok else "",
                             f"{s['E2']['velocity']:.3e}", f"{s['energyDrift']:+.2e}"])
                prev = e if ok else None
    header = ["scheme", "distortion", "cells", "nonOrthMax", "nonOrthAvg", "E2p",
              "q", "E2U", "energyDrift"]
    with (root / "studyWallLayer.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    widths = [max(len(str(v)) for v in col) for col in zip(header, *rows)]
    for row in [header, *rows]:
        print("  ".join(str(v).rjust(w) for v, w in zip(row, widths)))


if __name__ == "__main__":
    main()
