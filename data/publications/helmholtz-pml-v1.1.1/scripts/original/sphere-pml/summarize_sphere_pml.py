#!/usr/bin/env python3
"""Tables of the penetrable-sphere study.

Writes convergence.csv in each study directory; for the eight-rank study,
rank_comparison.csv with the relative L2 difference of the complex pressure
to the serial run of the same mesh in the main study; and for the main study,
rigid_difference.csv with the volume-weighted relative L2 difference between
the exact fields of the penetrable and a rigid sphere in the host cells of
the physical region on the finest exact-data mesh.
"""
import csv, json, sys
from pathlib import Path
import numpy as np
import sphere_reference
from sphere_pml import A, D, K, MATERIALS, field, physical

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
    main_rows = json.loads((root/"main/results.json").read_text()) if (root/"main/results.json").exists() else []
    finest = {}
    for r in main_rows:
        if r["completed"] and r["truncation"] == "dirichlet":
            finest[r["material"]] = max(finest.get(r["material"], 0), r["N"])
    if finest:
        out = []
        for material, n in finest.items():
            prepared = root/"main"/material/f"prepared-dirichlet-N{n}"
            C = field(prepared/"0/C"); V = field(prepared/"0/V", len(C))
            host = physical(C) & (np.linalg.norm(C, axis=1) > A)
            (rhog, cg), (rhol, cl) = MATERIALS[material]
            args = (C[host], A, D, K, K*cg/cl, rhog, rhol)
            p, rigid = sphere_reference.pressure(*args), sphere_reference.pressure(*args, rigid=True)
            out.append(dict(material=material, N=n, hostCells=int(host.sum()),
                            relL2Difference=float(np.sqrt(np.sum(V[host]*abs(p-rigid)**2)/np.sum(V[host]*abs(p)**2)))))
        with open(root/"main"/"rigid_difference.csv", "w", newline="") as f:
            w = csv.DictWriter(f, ["material", "N", "hostCells", "relL2Difference"])
            w.writeheader(); w.writerows(out)
        for r in out:
            print(f"{r['material']:<11} penetrable vs rigid sphere, host field: {r['relL2Difference']:.3e}")
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
