#!/usr/bin/env python3
"""Tabulate the two-layer cavity studies written by runStudy."""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np


def summary(case):
    return json.loads((case / "postProcessing/twoLayerVerification/"
                       "summary.json").read_text())


def probes(case):
    rows = [line.split() for line in
            (case / "postProcessing/cellProbes/0/p").read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    data = np.array(rows, dtype=float)
    return data[:, 0], data[:, 1:]


def fmt(x, digits=4):
    return "" if x is None else f"{x:.{digits}e}"


def order(coarse, fine):
    if not coarse or not fine or coarse <= 0 or fine <= 0:
        return ""
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


def study_a(root, out):
    rows = []
    for case in sorted((root / "A").iterdir()):
        s = summary(case)
        r = s["rates"]
        rows.append([
            case.name, s["setup"]["method"], s["setup"]["blockSolver"],
            f"{s['setup']['f']:.3f}",
            fmt(s["implementation"]["maxModalErrorVsDiscrete"], 2),
            fmt(r["gammaFit"]), fmt(r["gammaDiscrete"]),
            fmt(s["energy"]["discreteWaveEnergyMaxRelativeDrift"], 2),
        ])
    write(out / "studyA.csv",
          ["case", "method", "solver", "f", "maxModalErrDisc", "gammaFit",
           "gammaDisc", "EhDrift"], rows)


def refinement(root, out, letter):
    groups = {}
    for case in (root / letter).iterdir():
        s = summary(case)
        key = (s["setup"]["interfaceOffset"], s["setup"]["mode"])
        groups.setdefault(key, []).append((s, case))
    for (offset, mode), items in sorted(groups.items()):
        items.sort(key=lambda item: item[0]["setup"]["cellsPerUnit"])
        disc = [abs(s["rates"]["relativeFrequencyErrorDiscreteVsExact"])
                for s, _ in items]
        fit = [abs(s["rates"]["relativeFrequencyErrorFitVsExact"])
               for s, _ in items]
        ep = [s["E2"]["pressureVsExact"] for s, _ in items]
        eu = [s["E2"]["velocityVsExact"] for s, _ in items]
        rows = []
        for i, (s, _) in enumerate(items):
            prev = i - 1
            rows.append([
                f"{s['setup']['gasCellsPerWavelength']:.0f}", s["setup"]["nx"],
                fmt(s["rates"]["relativeFrequencyErrorDiscreteVsExact"]),
                fmt(fit[i]), order(disc[prev], disc[i]) if i else "",
                fmt(ep[i]), order(ep[prev], ep[i]) if i else "",
                fmt(eu[i]), order(eu[prev], eu[i]) if i else "",
                f"{s['initialModalShare']:.6f}",
                fmt(s["implementation"]["maxModalErrorVsDiscrete"], 1),
            ])
        write(out / f"study{letter}_s{offset:g}_{mode}.csv",
              ["N_lambda_g", "nx", "omegaH/omega-1", "|fit|", "order",
               "E2p", "order", "E2U", "order", "modalShare", "implErr"],
              rows)


def study_d(root, out):
    rows = []
    cases = {c.name: c for c in (root / "D").iterdir()}
    for name in sorted(n for n in cases if n.endswith("_np2")):
        serial, parallel = cases[name[:-4] + "_np1"], cases[name]
        s1, s2 = summary(serial), summary(parallel)
        t1, p1 = probes(serial)
        t2, p2 = probes(parallel)
        if not np.allclose(t1, t2):
            raise RuntimeError(f"probe times differ for {name}")
        proc = s2["processorFace"]
        rows.append([
            name[:-4], f"{proc['invRhofRank0']:.6g}",
            f"{proc['invRhofRank1']:.6g}", s2["operatorSymmetric"],
            fmt(s1["rates"]["relativeFrequencyErrorDiscreteVsExact"]),
            fmt(s2["rates"]["relativeFrequencyErrorDiscreteVsExact"]),
            fmt(s2["energy"]["discreteWaveEnergyMaxRelativeDrift"], 2),
            fmt(float(np.abs(p1 - p2).max()), 2),
        ])
    write(out / "studyD.csv",
          ["case", "invRhofRank0", "invRhofRank1", "symmetric",
           "freqErrSerial", "freqErrParallel", "EhDriftParallel",
           "maxDiffSerialParallel"], rows)


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "studies")
    if (root / "A").is_dir():
        study_a(root, root)
    for letter in ("B", "C"):
        if (root / letter).is_dir():
            refinement(root, root, letter)
    if (root / "D").is_dir():
        study_d(root, root)


if __name__ == "__main__":
    main()
