#!/usr/bin/env python3
"""Tabulate the cost study written by runCostStudy.

For each run: solver CPU time (ExecutionTime at the end), set-up time
(ExecutionTime after the first step, including assembly and factorization
of the monolithic matrix), time per step, E2 of the pressure against the
free-space solution, and, on h = 1 mm, the temporal E2 against the
reference run at Co = 1/16 on the same mesh, over the common samples.
"""
import csv
import json
import math
import re
import sys
from pathlib import Path

import numpy as np


def execution_times(log):
    return [float(x) for x in re.findall(r"ExecutionTime = ([0-9.eE+-]+) s", log.read_text())]


def pressure(case):
    rows = [line.split() for line in
            (case / "postProcessing/cellProbes/0/p").read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    data = np.array(rows, dtype=float)
    return data[:, 0], data[:, 1:]


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "studies/cost")
    rows = []
    for mesh in sorted(d for d in root.iterdir() if d.is_dir()):
        ref = mesh / "reference_Co0.0625"
        t_ref, p_ref = pressure(ref) if ref.is_dir() else (None, None)
        for case in sorted(mesh.iterdir()):
            log = case / "log.acousticWaveFoam"
            path = case / "postProcessing/pmlVerification/summary.json"
            timing = case.name.startswith("timing_")
            if not log.is_file() or not (timing or path.is_file()):
                continue
            s = json.loads((case / "caseSetup.json").read_text())
            times = execution_times(log)
            summ = {} if timing else json.loads(path.read_text())
            temporal = ""
            if t_ref is not None and case != ref and not timing:
                t, p = pressure(case)
                # Common samples by index: the sample interval is 4 h/c.
                period = 4 * s["h"] / s["cG"]
                index = np.rint(t / period).astype(int)
                index_ref = np.rint(t_ref / period).astype(int)
                common = np.intersect1d(index, index_ref)
                keep = np.isin(index, common)
                keep_ref = np.isin(index_ref, common)
                if keep.sum() and np.isfinite(p[keep]).all():
                    temporal = f"{math.sqrt(np.sum((p[keep] - p_ref[keep_ref])**2) / np.sum(p_ref[keep_ref]**2)):.3e}"
                elif keep.sum():
                    temporal = "unstable"
            update = case.name.split("_Co")[0]
            setup = times[0] if times else float("nan")
            total = times[-1] if times else float("nan")
            steps = s["nSteps"]
            rows.append([mesh.name[1:], update, f"{s['courant']:g}", steps,
                         f"{setup:.3g}", f"{total:.3g}",
                         f"{(total - setup) / max(steps - 1, 1) * 1e3:.3g}",
                         f"{summ['E2exact']['pressure']:.3e}" if summ.get("E2exact") else "",
                         temporal, summ.get("finite", "")])
    header = ["h", "update", "Co", "steps", "setup_s", "total_s", "ms_per_step",
              "E2p_exact", "E2p_temporal", "finite"]
    with (root / "studyCost.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    widths = [max(len(str(v)) for v in col) for col in zip(header, *rows)]
    for row in [header, *rows]:
        print("  ".join(str(v).rjust(w) for v, w in zip(row, widths)))


if __name__ == "__main__":
    main()
