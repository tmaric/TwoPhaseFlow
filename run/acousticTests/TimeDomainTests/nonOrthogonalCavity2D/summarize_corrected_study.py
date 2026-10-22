#!/usr/bin/env python3
"""Tabulate the corrected-Laplacian study written by runCorrectedStudy.

S: E2 against the exact mode and the energy drift over 10 periods.
T: temporal E2 of each run against the run with STEPS_PER_PERIOD 2560 of
the same update on the same mesh (common samples, 8 per period).
L: long-time growth and energy drift over 200 periods.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np


def summary(case):
    path = case / "postProcessing/cavity2dVerification/summary.json"
    return json.loads(path.read_text()) if path.is_file() else None


def probes(case):
    rows = [line.split() for line in
            (case / "postProcessing/cellProbes/0/p").read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    data = np.array(rows, dtype=float)
    return data[:, 0], data[:, 1:]


def write(path, header, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"\n{path}")
    widths = [max(len(str(v)) for v in col) for col in zip(header, *rows)]
    for row in [header, *rows]:
        print("  ".join(str(v).rjust(w) for v, w in zip(row, widths)))


def fmt(x):
    return "" if x is None or not math.isfinite(x) else f"{x:.3e}"


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "studies/corrected")
    if (root / "S").is_dir():
        rows = []
        for update in ("monolithic", "segregated"):
            for eps in ("0.15", "0.3", "0.6"):
                prev = None
                for n in (20, 40, 80, 160):
                    s = summary(root / "S" / f"{update}_eps{eps}_N{n}")
                    if s is None:
                        continue
                    e = s["E2"]["pressure"]
                    ok = math.isfinite(e) and e < 10
                    q = (f"{math.log(prev / e, 2):.2f}" if prev and ok else "")
                    rows.append([update, eps, n, fmt(e), q, fmt(s["E2"]["velocity"]),
                                 f"{s['energyDrift']:+.2e}"])
                    prev = e if ok else None
        write(root / "studyS.csv", ["update", "distortion", "cells", "E2p", "q",
                                    "E2U", "energyDrift"], rows)
    if (root / "T").is_dir():
        rows = []
        for update in ("monolithic", "segregated"):
            ref = root / "T" / f"{update}_NT2560"
            t_ref, p_ref = probes(ref)
            prev = None
            for nt in (40, 80, 160, 320, 640):
                case = root / "T" / f"{update}_NT{nt}"
                s = summary(case)
                t, p = probes(case)
                period = 1 / s["setup"]["f"] / 8
                i, i_ref = np.rint(t / period).astype(int), np.rint(t_ref / period).astype(int)
                common = np.intersect1d(i, i_ref)
                a, b = p[np.isin(i, common)], p_ref[np.isin(i_ref, common)]
                e = math.sqrt(np.sum((a - b) ** 2) / np.sum(b**2))
                rows.append([update, nt, f"{s['setup']['courant']:.3f}", fmt(e),
                             f"{math.log(prev / e, 2):.2f}" if prev else "",
                             fmt(s["E2"]["pressure"]), f"{s['energyDrift']:+.2e}"])
                prev = e
        write(root / "studyT.csv", ["update", "N_T", "Co", "E2p_temporal", "q",
                                    "E2p_exact", "energyDrift"], rows)
    if (root / "L").is_dir():
        rows = []
        for case in sorted((root / "L").iterdir()):
            s = summary(case)
            if s is None:
                rows.append([case.name, "", "", "no summary", "", ""])
                continue
            rows.append([case.name, f"{s['setup']['courant']:.2f}",
                         f"{s['lastTime']:.4g}", s["finite"],
                         f"{s['modal']['growthPerPeriod'] - 1:+.3e}",
                         f"{s['energyDrift']:+.3e}"])
        write(root / "studyL.csv", ["case", "Co", "tLast", "finite",
                                    "modalGrowthPerPeriod", "energyDrift"], rows)
    if (root / "N").is_dir():
        corrector_scan(root)


def corrector_scan(root):
    import re
    rows = []
    for case in sorted((root / "N").iterdir(), key=lambda c: int(c.name.split("corr")[1])):
        s = summary(case)
        log = (case / "log.acousticWaveFoam").read_text()
        changes = re.findall(r"Non-orthogonal corrector (\d+): relative pressure change = ([0-9.eE+-]+)", log)
        last = {}
        for k, v in changes:
            last[int(k)] = float(v)
        times = [float(x) for x in re.findall(r"ExecutionTime = ([0-9.eE+-]+) s", log)]
        steps = s["setup"]["nSteps"] if s else 1
        rows.append([case.name.split("corr")[1],
                     fmt(s["E2"]["pressure"]) if s else "",
                     f"{s['energyDrift']:+.2e}" if s else "",
                     fmt(last[max(last)]) if last else "",
                     f"{(times[-1] - times[0]) / max(steps - 1, 1) * 1e3:.3g}" if times else ""])
    write(root / "studyN.csv", ["correctors", "E2p", "energyDrift",
                                "lastCorrectorChange", "ms_per_step"], rows)


if __name__ == "__main__":
    main()
