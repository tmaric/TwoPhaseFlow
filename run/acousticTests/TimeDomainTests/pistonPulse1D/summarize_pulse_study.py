#!/usr/bin/env python3
"""Tabulate the piston-pulse studies written by runStudy."""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np


def pressure_history(case):
    s = json.loads((case / "caseSetup.json").read_text())
    rows = [line.split() for line in
            (case / "postProcessing/cellProbes/0/p").read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    data = np.array(rows, dtype=float)
    t, p = data[:, 0], data[:, 1:]
    if t[0] > 0.5 * s["deltaT"]:
        t = np.concatenate(([0.0], t))
        p = np.vstack((np.zeros(s["nx"]), p))
    return s, t, p


def temporal_errors(root, cases):
    """E2 of each run against the run with N_T = 1400 on the same mesh.

    The reference's temporal error is negligible, so the difference isolates
    the temporal error of the coarser time step.
    """
    reference = root / "B" / "NX2800"
    if not reference.is_dir():
        return {}
    s, t_ref, p_ref = pressure_history(reference)
    x = (np.arange(s["nx"]) + 0.5) * s["h"]
    beta = 1.0 / np.where(x < s["xInterface"], s["rho1"] * s["c1"] ** 2,
                          s["rho2"] * s["c2"] ** 2)
    weight = np.where(x < s["pmlStart"], beta, 0.0)[None, :]
    errors = {}
    for case in cases:
        _, t, p = pressure_history(case)
        if not np.allclose(t, t_ref):
            raise RuntimeError(f"sample times differ for {case}")
        errors[case.name] = math.sqrt(np.sum(weight * (p - p_ref) ** 2)
                                      / np.sum(weight * p_ref**2))
    return errors


def summary(case):
    return json.loads((case / "postProcessing/pulseVerification/"
                       "summary.json").read_text())


def fmt(x, digits=3):
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


def convergence(items, key):
    rows = []
    for i, s in enumerate(items):
        e, prev = s["E2"], items[i - 1]["E2"] if i else None
        rows.append([
            key(s), f"{s['setup']['gasCellsPerWavelength']:.1f}",
            s["setup"]["stepsPerPeriod"],
            fmt(e["pressure"]), order(prev["pressure"], e["pressure"]) if i else "",
            fmt(e["velocity"]), order(prev["velocity"], e["velocity"]) if i else "",
            fmt(e["pressureMedium2"]),
            order(prev["pressureMedium2"], e["pressureMedium2"]) if i else "",
            fmt(s["pml"]["spuriousRelativePressure"]),
        ])
    return rows


HEADER = ["case", "N_lambda_g", "N_T", "E2p", "order", "E2U", "order",
          "E2p_medium2", "order", "spurious"]


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "studies")
    if (root / "A").is_dir():
        for method in ("trapezoidal", "backward"):
            cases = sorted((root / "A").glob(f"{method}_*"),
                           key=lambda c: summary(c)["setup"]["stepsPerPeriod"])
            if not cases:
                continue
            temporal = temporal_errors(root, cases)
            rows = convergence([summary(c) for c in cases],
                               lambda s: s["setup"]["method"])
            for i, (row, case) in enumerate(zip(rows, cases)):
                e = temporal.get(case.name)
                prev = temporal.get(cases[i - 1].name) if i else None
                row += [fmt(e), order(prev, e) if prev else ""]
            write(root / f"studyA_{method}.csv",
                  HEADER + ["E2p_temporal", "order"], rows)
    for letter in ("B", "D"):
        if (root / letter).is_dir():
            items = sorted((summary(c) for c in (root / letter).iterdir()),
                           key=lambda s: s["setup"]["nx"])
            write(root / f"study{letter}.csv", HEADER,
                  convergence(items, lambda s: f"{s['setup']['config']} "
                                               f"NX{s['setup']['nx']}"))
    if (root / "C").is_dir():
        items = sorted((summary(c) for c in (root / "C").iterdir()),
                       key=lambda s: (s["setup"]["pmlThickness"],
                                      s["setup"]["sigmaMax"]))
        rows = [[s["setup"]["pmlThickness"], f"{s['setup']['sigmaMax']:g}",
                 fmt(math.exp(-2 * s["setup"]["sigmaMax"]
                              * s["setup"]["pmlThickness"]
                              / ((s["setup"]["pmlOrder"] + 1) * s["setup"]["c2"]))),
                 fmt(s["pml"]["spuriousRelativePressure"]),
                 fmt(s["E2"]["pressureMedium2"]), fmt(s["E2"]["pressure"])]
                for s in items]
        write(root / "studyC.csv",
              ["thickness", "sigmaMax", "continuousRoundTrip", "spurious",
               "E2p_medium2", "E2p"], rows)


if __name__ == "__main__":
    main()
