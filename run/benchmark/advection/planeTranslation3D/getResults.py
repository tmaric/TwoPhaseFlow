#!/usr/bin/env python3
"""Collect the plane-translation gate and check it.

Reads, per case under Cases/<velocity>/<davofMethod>/N<N>/, the per-step
areaFractionError.dat (L1, L2, Linf of |alphaf - alphaf_exact| over the faces
of the interface cells, the face and cell counts, the consistency
max |m_c - A n|/V^(2/3)), the per-step reconstructionError error.dat
(LNormalDiff1/2/Inf, LCentreDiff1/2/Inf against the translated plane) and
volumeFractionError.dat (E1 = sum |alpha - alpha_exact| V), and writes

    results_planeTranslation3D.csv   one row per case: the maxima over the steps
    history_planeTranslation3D.csv   one row per case and step

Gate (every case): max over the steps of L1(alphaf) < GATE and of
LNormalDiff1 < GATE (1e-12); the maxima (Linf) are reported; the script exits
non-zero on a failure.
"""
import csv
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NAME = os.path.basename(HERE)
GATE = 1e-12


def read_param(path):
    vals = {}
    for line in open(path):
        m = re.match(r"^\s*([A-Za-z_]\w*)\s+([^;]+);", line)
        if m:
            vals[m.group(1)] = m.group(2).strip()
    return vals


def read_dat(path):
    """All rows of a writeFile table as dicts keyed by the header names."""
    names, rows = None, []
    if not os.path.isfile(path):
        return rows
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            if s[1:].split() and s[1:].split()[0] == "Time":
                names = s[1:].split()
            continue
        if names:
            rows.append(dict(zip(names, s.split())))
    return rows


def main():
    results, history = [], []
    for case in sorted(glob.glob(os.path.join(HERE, "Cases", "*", "*", "N*"))):
        par = read_param(os.path.join(case, "system", "simulationParameter"))
        af = read_dat(os.path.join(case, "areaFractionError.dat"))
        rec = read_dat(os.path.join(case, "error.dat"))
        vol = read_dat(os.path.join(case, "volumeFractionError.dat"))
        if not af:
            print(f"[getResults] no areaFractionError.dat in {case}")
            continue
        n, L = int(par["N"]), float(par["L"])
        h = L/n
        rec_by_t = {r["Time"]: r for r in rec}
        vol_by_t = {r["Time"]: r for r in vol}
        label = {"velocity": par["VELOCITY"], "davofMethod": par["DAVOFMETHOD"], "N": n, "h": h,
                 "case": os.path.relpath(case, HERE)}
        maxima = {"L1": 0.0, "L2": 0.0, "Linf": 0.0, "consistency": 0.0,
                  "LNormalDiff1": 0.0, "LNormalDiffInf": 0.0, "LCentreDiff2_over_h": 0.0,
                  "LCentreDiffInf_over_h": 0.0, "E1": 0.0}
        for r in af:
            t = r["Time"]
            rr = rec_by_t.get(t, {})
            vv = vol_by_t.get(t, {})
            row = dict(label)
            row.update({"Time": t, "L1": float(r["L1"]), "L2": float(r["L2"]), "Linf": float(r["Linf"]),
                        "nFaces": int(r["nFaces"]), "nInterfaceCells": int(r["nInterfaceCells"]),
                        "consistency": float(r["consistency"]),
                        "LNormalDiff1": float(rr.get("LNormalDiff1", "nan")),
                        "LNormalDiffInf": float(rr.get("LNormalDiffInf", "nan")),
                        "LCentreDiff2_over_h": float(rr.get("LCentreDiff2", "nan"))/h,
                        "LCentreDiffInf_over_h": float(rr.get("LCentreDiffInf", "nan"))/h,
                        "E1": float(vv.get("Eshape", "nan"))})
            history.append(row)
            for k in maxima:
                v = row[k]
                if v == v:  # not nan
                    maxima[k] = max(maxima[k], v)
        res = dict(label)
        res["nSteps"] = len(af)
        res.update({f"max_{k}": v for k, v in maxima.items()})
        res["nInterfaceCells_last"] = int(af[-1]["nInterfaceCells"])
        res["PASS"] = maxima["L1"] < GATE and maxima["LNormalDiff1"] < GATE
        results.append(res)
    if not results:
        print("[getResults] nothing found under Cases/")
        return 1

    for name, rows in (("results", results), ("history", history)):
        path = os.path.join(HERE, f"{name}_{NAME}.csv")
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"[getResults] wrote {path} ({len(rows)} rows)")

    status = 0
    print(f"[getResults] gate (GATE = {GATE:g} on max over the steps of L1(alphaf) and L1(n)):")
    print(f"[getResults] {'velocity':8s} {'method':15s} N    L1(af)   Linf(af)  L1(n)    Linf(n)   L2(x)/h   Linf(x)/h  consist   E1        cells  steps")
    for r in results:
        status |= (0 if r["PASS"] else 2)
        print(f"[getResults] {r['velocity']:8s} {r['davofMethod']:15s} {r['N']:<4} {r['max_L1']:.2e} {r['max_Linf']:.2e}  "
              f"{r['max_LNormalDiff1']:.2e} {r['max_LNormalDiffInf']:.2e}  {r['max_LCentreDiff2_over_h']:.2e}  "
              f"{r['max_LCentreDiffInf_over_h']:.2e}   {r['max_consistency']:.2e}  {r['max_E1']:.2e}  "
              f"{r['nInterfaceCells_last']:<5} {r['nSteps']:<3} {'PASS' if r['PASS'] else 'FAIL'}")
    return status


if __name__ == "__main__":
    sys.exit(main())
