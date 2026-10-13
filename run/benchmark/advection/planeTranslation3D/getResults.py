#!/usr/bin/env python3
"""Collect the plane-translation gate and check it.

Reads, per case under Cases/<velocity>/<row>/N<N>/, the per-step
areaFractionError.dat (L1, L2, Linf of |alphaf - alphaf_exact| over the faces
of the interface cells, the face and cell counts, the consistency
max |m_c - A n|/V^(2/3); -1 for the comparator rows without face fractions),
the per-step reconstructionError error.dat (LNormalDiff1/2/Inf,
LCentreDiff1/2/Inf against the translated plane), volumeFractionError.dat
(E1 = sum |alpha - alpha_exact| V) and the per-step sliver count of
log.advectorVoF, and writes

    results_planeTranslation3D.csv   one row per case: the row parameters
                                     (scheme, davofMethod, sliverNormal,
                                     sliverTol, normalSource), the maxima
                                     over the steps and the last step
    history_planeTranslation3D.csv   one row per case and step
    growth_planeTranslation3D.csv    one row per case: L1(alphaf) of the
                                     first 20 steps (L1_1 ... L1_20), the
                                     per-step growth factors
                                     g_k = L1(k)/L1(k-1) (g_2 ... g_20),
                                     their maximum and geometric mean

Gate (the davof rows): max over the steps of L1(alphaf) < GATE and of
LNormalDiff1 < GATE (1e-12); the maxima (Linf) are reported; the script
exits non-zero on a failure. The comparator rows are reported, not gated.
"""
import csv
import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NAME = os.path.basename(HERE)
GATE = 1e-12
NHIST = 20
PARAMS = ["RECONSCHEME", "DAVOFMETHOD", "SLIVERNORMAL", "SLIVERTOL", "NORMALSOURCE"]


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


def read_slivers(path):
    """Per time step the (interface cells, slivers) of the advection step's
    own reconstruction: the first davof Info line after each 'Time = ' line
    (the function objects reconstruct again at the end of the step)."""
    out = []
    if not os.path.isfile(path):
        return out
    expect = False
    for line in open(path):
        if line.startswith("Time = "):
            expect = True
            continue
        m = re.match(r"^davof: (\d+) interface cells .*?, (\d+) slivers", line)
        if m and expect:
            out.append((int(m.group(1)), int(m.group(2))))
            expect = False
    return out


def growth_row(label, l1):
    """L1_1..L1_NHIST, g_2..g_NHIST, max and geometric mean of the factors."""
    row = dict(label)
    l1 = [v if (v == v and v >= 0) else float("nan") for v in l1[:NHIST]]
    for k in range(1, NHIST + 1):
        row[f"L1_{k}"] = l1[k - 1] if k <= len(l1) else float("nan")
    gs = []
    for k in range(2, NHIST + 1):
        g = float("nan")
        if k <= len(l1) and l1[k - 2] > 0 and l1[k - 1] == l1[k - 1]:
            g = l1[k - 1]/l1[k - 2]
        row[f"g_{k}"] = g
        if g == g and g > 0:
            gs.append(g)
    row["maxGrowth"] = max(gs) if gs else float("nan")
    row["geoMeanGrowth"] = math.exp(sum(math.log(g) for g in gs)/len(gs)) if gs else float("nan")
    return row


def main():
    results, history, growth = [], [], []
    for case in sorted(glob.glob(os.path.join(HERE, "Cases", "*", "*", "N*"))):
        par = read_param(os.path.join(case, "system", "simulationParameter"))
        af = read_dat(os.path.join(case, "areaFractionError.dat"))
        rec = read_dat(os.path.join(case, "error.dat"))
        vol = read_dat(os.path.join(case, "volumeFractionError.dat"))
        slv = read_slivers(os.path.join(case, "log.advectorVoF"))
        if not af:
            print(f"[getResults] no areaFractionError.dat in {case}")
            continue
        n, L = int(par["N"]), float(par["L"])
        h = L/n
        davof = par["RECONSCHEME"] == "davof"
        rec_by_t = {r["Time"]: r for r in rec}
        vol_by_t = {r["Time"]: r for r in vol}
        label = {"velocity": par["VELOCITY"], "row": os.path.basename(os.path.dirname(case)),
                 "scheme": par["RECONSCHEME"], "davofMethod": par.get("DAVOFMETHOD", ""),
                 "sliverNormal": par.get("SLIVERNORMAL", ""), "sliverTol": par.get("SLIVERTOL", ""),
                 "normalSource": par.get("NORMALSOURCE", ""), "N": n, "h": h,
                 "case": os.path.relpath(case, HERE)}
        maxima = {"L1": 0.0, "L2": 0.0, "Linf": 0.0, "consistency": 0.0,
                  "LNormalDiff1": 0.0, "LNormalDiffInf": 0.0, "LCentreDiff2_over_h": 0.0,
                  "LCentreDiffInf_over_h": 0.0, "E1": 0.0}
        l1hist = []
        for k, r in enumerate(af):
            t = r["Time"]
            rr = rec_by_t.get(t, {})
            vv = vol_by_t.get(t, {})
            row = dict(label)
            # the sliver count of the advection step's own reconstruction
            ns = slv[k] if k < len(slv) else (None, None)
            row.update({"step": k + 1, "Time": t, "L1": float(r["L1"]), "L2": float(r["L2"]),
                        "Linf": float(r["Linf"]), "nFaces": int(r["nFaces"]),
                        "nInterfaceCells": int(r["nInterfaceCells"]),
                        "consistency": float(r["consistency"]),
                        "nSliver": ns[1] if ns[1] is not None else "",
                        "LNormalDiff1": float(rr.get("LNormalDiff1", "nan")),
                        "LNormalDiffInf": float(rr.get("LNormalDiffInf", "nan")),
                        "LCentreDiff2_over_h": float(rr.get("LCentreDiff2", "nan"))/h,
                        "LCentreDiffInf_over_h": float(rr.get("LCentreDiffInf", "nan"))/h,
                        "E1": float(vv.get("Eshape", "nan"))})
            history.append(row)
            l1hist.append(row["L1"])
            for key in maxima:
                v = row[key]
                if v == v:  # not nan
                    maxima[key] = max(maxima[key], v)
        res = dict(label)
        res["nSteps"] = len(af)
        res.update({f"max_{k}": v for k, v in maxima.items()})
        res["L1_last"] = l1hist[-1]
        res["nInterfaceCells_last"] = int(af[-1]["nInterfaceCells"])
        res["nSliver_max"] = max((s[1] for s in slv), default="")
        res["PASS"] = (maxima["L1"] < GATE and maxima["LNormalDiff1"] < GATE) if davof else ""
        results.append(res)
        if davof:
            growth.append(growth_row(label, l1hist))
    if not results:
        print("[getResults] nothing found under Cases/")
        return 1

    for name, rows in (("results", results), ("history", history), ("growth", growth)):
        if not rows:
            continue
        path = os.path.join(HERE, f"{name}_{NAME}.csv")
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"[getResults] wrote {path} ({len(rows)} rows)")

    status = 0
    print(f"[getResults] gate (GATE = {GATE:g} on max over the steps of L1(alphaf) and L1(n), davof rows):")
    print(f"[getResults] {'velocity':8s} {'row':36s} N    L1(af)   Linf(af)  L1(n)    Linf(n)   L2(x)/h   Linf(x)/h  consist   E1        cells  slv  steps")
    for r in results:
        if r["PASS"] == "":
            verdict = "-"
        else:
            verdict = "PASS" if r["PASS"] else "FAIL"
            status |= (0 if r["PASS"] else 2)
        print(f"[getResults] {r['velocity']:8s} {r['row']:36s} {r['N']:<4} {r['max_L1']:.2e} {r['max_Linf']:.2e}  "
              f"{r['max_LNormalDiff1']:.2e} {r['max_LNormalDiffInf']:.2e}  {r['max_LCentreDiff2_over_h']:.2e}  "
              f"{r['max_LCentreDiffInf_over_h']:.2e}   {r['max_consistency']:.2e}  {r['max_E1']:.2e}  "
              f"{r['nInterfaceCells_last']:<5} {str(r['nSliver_max']):4s} {r['nSteps']:<3} {verdict}")
    if growth:
        print(f"[getResults] L1(alphaf) history (first {NHIST} steps) and the growth factors L1(k)/L1(k-1), davof rows:")
        for g in growth:
            l1 = [g[f"L1_{k}"] for k in range(1, NHIST + 1) if g[f"L1_{k}"] == g[f"L1_{k}"]]
            gs = [g[f"g_{k}"] for k in range(2, NHIST + 1) if g[f"g_{k}"] == g[f"g_{k}"]]
            print(f"[getResults] {g['velocity']} {g['row']} N{g['N']}: max growth {g['maxGrowth']:.2f}, "
                  f"geometric mean {g['geoMeanGrowth']:.2f}")
            print("[getResults]   L1: " + " ".join(f"{v:.1e}" for v in l1))
            print("[getResults]   g:  " + " ".join(f"{v:.2f}" for v in gs))
    return status


if __name__ == "__main__":
    sys.exit(main())
