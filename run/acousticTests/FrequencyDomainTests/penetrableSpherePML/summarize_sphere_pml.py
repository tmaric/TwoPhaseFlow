#!/usr/bin/env python3
"""Tables of the penetrable-sphere study.

Writes convergence.csv in each study directory and, for the eight-rank
study, rank_comparison.csv with the relative L2 difference of the complex
pressure to the serial run of the same mesh in the main study.
"""
import csv, json, sys
from pathlib import Path
import numpy as np
from sphere_pml import field

COLUMNS = ["material", "truncation", "N", "cells", "physicalCells", "physicalRelL2", "observedOrder",
           "hostRelL2", "dropRelL2", "interfaceRelL2", "farHostRelL2", "minusDirichletRelL2",
           "sigmaOverOmega", "pmlLength", "ranks", "seconds", "maxRelativeResidual"]

def pressure(case):
    return field(case/"1/Pre") + 1j*field(case/"1/Pim")

def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "studyResults")
    failures = 0
    for study in sorted(p.parent for p in root.glob("*/results.json")):
        rows = json.loads((study/"results.json").read_text())
        failures += sum(not r["completed"] for r in rows)
        with open(study/"convergence.csv", "w", newline="") as f:
            w = csv.DictWriter(f, COLUMNS, extrasaction="ignore")
            w.writeheader(); w.writerows(rows)
        print(f"== {study.name}")
        for r in rows:
            if r["completed"]:
                order = f"{r['observedOrder']:.2f}" if "observedOrder" in r else ""
                print(f"{r['material']:<11} {r['truncation']:<13} N={r['N']:<3} "
                      f"{100*r['physicalRelL2']:.4f} % {order}")
    ranks8 = root/"ranks8"
    if (ranks8/"results.json").exists():
        out = []
        for r in json.loads((ranks8/"results.json").read_text()):
            serial = root/"main"/r["material"]/f"{r['truncation']}-N{r['N']}"
            if not r["completed"] or not serial.exists():
                continue
            p1, p8 = pressure(serial), pressure(ranks8/r["material"]/f"{r['truncation']}-N{r['N']}")
            out.append(dict(material=r["material"], truncation=r["truncation"], N=r["N"], ranks=r["ranks"],
                            relL2Difference=float(np.linalg.norm(p8-p1)/np.linalg.norm(p1))))
        with open(ranks8/"rank_comparison.csv", "w", newline="") as f:
            w = csv.DictWriter(f, ["material", "truncation", "N", "ranks", "relL2Difference"])
            w.writeheader(); w.writerows(out)
        for r in out:
            print(f"{r['material']:<11} {r['truncation']:<13} N={r['N']:<3} "
                  f"{r['ranks']} ranks vs serial: {r['relL2Difference']:.2e}")
    if failures:
        raise SystemExit(f"{failures} failed runs")

if __name__ == "__main__":
    main()
