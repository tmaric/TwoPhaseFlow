#!/usr/bin/env python3
"""Tabulate the amplification-matrix study written by runSpectralStudy."""
import csv
import json
import sys
from pathlib import Path


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "studies/C")
    rows = []
    for case in sorted(root.iterdir()):
        path = case / "spectrum.json"
        if not path.is_file():
            continue
        r = json.loads(path.read_text())
        update = ("monolithic" if r["method"] == "blockCoupled" else
                  "segregated, implicit damping" if r["segregatedImplicitDamping"]
                  else "segregated")
        pn = r["powerNorm"]
        rows.append([
            case.name, update, f"{r['courant']:g}", f"{r['sigmaMax']:.4g}",
            f"{r['outerCellSigmaDt']:.3f}", f"{r['spectralRadius'] - 1:.3e}",
            r["nUnstable"], f"{r['explicitDampingRadius'] - 1:.3e}",
            f"{r['nearOne']['algebraic']}/{r['nearOne']['geometric']}",
            f"{pn['last']:.3e}", pn["lastN"],
            f"{r['superpositionError']:.1e}",
            f"{r.get('cavityEigenvalueError', float('nan')):.1e}",
        ])
    header = ["case", "update", "Co", "sigmaMax", "outerSigmaDt", "rho-1",
              "nUnstable", "explicitRadius-1", "g=1 alg/geo", "normLast",
              "nLast", "superposition", "cavityEig"]
    with (root / "studyC.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    widths = [max(len(str(v)) for v in col) for col in zip(header, *rows)]
    for row in [header, *rows]:
        print("  ".join(str(v).rjust(w) for v, w in zip(row, widths)))


if __name__ == "__main__":
    main()
