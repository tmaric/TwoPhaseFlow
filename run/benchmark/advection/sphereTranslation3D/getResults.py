#!/usr/bin/env python3
"""Collect the sphere-translation ladder into four CSVs.

Per case under Cases/<row>/N<N>/: the last row of areaFractionError.dat
(L1, L2, Linf of |alphaf - alphaf_exact| over the faces of the interface
cells, -1 for the rows without face fractions; the consistency
max |m_c - A n|/V^(2/3)), of error.dat (the reconstruction error at the end
time: LNormalDiff1/2/Inf, LCentreDiff1/2/Inf against the translated sphere)
and of volumeFractionError.dat (Eshape = E1 = sum |alpha - alpha_exact| V,
Emass, Ebound), plus the last "t = ..." line of log.advectorVoF (sum(alpha V),
1 - max(alpha), min(alpha)) and the sliver counts of the davof Info lines.
Writes next to this script:

    results_sphereTranslation3D.csv   one row per case: the row parameters
                                      (scheme, davofMethod, sliverNormal,
                                      sliverTol, normalSource) and the
                                      end-of-run numbers
    orders_sphereTranslation3D.csv    pairwise orders over N per row for
                                      L1(alphaf), LNormalDiff1/2, LCentreDiff1/2, E1
    history_sphereTranslation3D.csv   the per-step areaFractionError rows of
                                      every davof case with the sliver count
    growth_sphereTranslation3D.csv    one row per davof case: L1(alphaf) of
                                      the first 20 steps (L1_1 ... L1_20), the
                                      per-step growth factors
                                      g_k = L1(k)/L1(k-1) (g_2 ... g_20),
                                      their maximum and geometric mean

The prediction of the proposal (finding (v), the iterated regeneration): the
quadratic initial state loses its second order within the first steps (the
plane clipping is q = 1), so L1(alphaf) after one diameter of travel is first
order, the normal first order, the centre second order.
"""
import csv
import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NAME = os.path.basename(HERE)
METRICS = ["L1", "L2", "LNormalDiff1", "LNormalDiff2", "LCentreDiff1", "LCentreDiff2", "E1"]
NHIST = 20


def read_param(path):
    vals = {}
    for line in open(path):
        m = re.match(r"^\s*([A-Za-z_]\w*)\s+([^;]+);", line)
        if m:
            vals[m.group(1)] = m.group(2).strip()
    return vals


def read_dat(path):
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


def read_log(path):
    """The last 't = ..., sum(alpha*V) = ..., 1-max(alpha1) = ..., min(alpha1) = ...' line."""
    out = {}
    if not os.path.isfile(path):
        return out
    for line in open(path):
        m = re.match(r"^t = ([^,]+),\s*sum\(alpha\*V\) = ([^,]+),\s*dev = ([^%]+)%,\s*1-max\(alpha1\) = (\S+) .*min\(alpha1\) = (\S+) ", line)
        if m:
            out = {"sumAlphaV": float(m.group(2)), "dev_percent": float(m.group(3)),
                   "oneMinusMaxAlpha": float(m.group(4)), "minAlpha": float(m.group(5))}
    return out


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


def fnum(d, k):
    try:
        return float(d[k])
    except (KeyError, ValueError, TypeError):
        return float("nan")


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
    for case in sorted(glob.glob(os.path.join(HERE, "Cases", "*", "N*"))):
        par = read_param(os.path.join(case, "system", "simulationParameter"))
        af = read_dat(os.path.join(case, "areaFractionError.dat"))
        rec = read_dat(os.path.join(case, "error.dat"))
        vol = read_dat(os.path.join(case, "volumeFractionError.dat"))
        log = read_log(os.path.join(case, "log.advectorVoF"))
        slv = read_slivers(os.path.join(case, "log.advectorVoF"))
        if not af and not rec:
            print(f"[getResults] no results in {case}")
            continue
        n = int(par["N"])
        row = os.path.basename(os.path.dirname(case))
        label = {"row": row, "scheme": par["RECONSCHEME"], "davofMethod": par.get("DAVOFMETHOD", ""),
                 "sliverNormal": par.get("SLIVERNORMAL", ""), "sliverTol": par.get("SLIVERTOL", ""),
                 "normalSource": par.get("NORMALSOURCE", ""), "N": n, "h": 1.0/n,
                 "case": os.path.relpath(case, HERE)}
        r = dict(label)
        r["nSteps"] = len(af)
        last = af[-1] if af else {}
        r.update({"L1": fnum(last, "L1"), "L2": fnum(last, "L2"), "Linf": fnum(last, "Linf"),
                  "nFaces": last.get("nFaces", ""), "nInterfaceCells": last.get("nInterfaceCells", ""),
                  "consistency": fnum(last, "consistency")})
        lr = rec[-1] if rec else {}
        for k in ("LNormalDiff1", "LNormalDiff2", "LNormalDiffInf", "LCentreDiff1", "LCentreDiff2",
                  "LCentreDiffInf", "nInterface", "areaSum", "areaExact"):
            r[k] = fnum(lr, k)
        lv = vol[-1] if vol else {}
        r["E1"] = fnum(lv, "Eshape")
        r["Emass"] = fnum(lv, "Emass")
        r["Ebound"] = fnum(lv, "Ebound")
        r.update(log)
        r["nSliver_max"] = max((s[1] for s in slv), default="")
        r["nSliver_last"] = slv[-1][1] if slv else ""
        results.append(r)
        if par["RECONSCHEME"] == "davof":
            l1hist = []
            for k, s in enumerate(af):
                ns = slv[k] if k < len(slv) else (None, None)
                history.append({"row": row, "N": n, "step": k + 1, **{kk: s[kk] for kk in s},
                                "nSliver": ns[1] if ns[1] is not None else ""})
                l1hist.append(fnum(s, "L1"))
            growth.append(growth_row(label, l1hist))
    if not results:
        print("[getResults] nothing found under Cases/")
        return 1

    def write(name, rows):
        path = os.path.join(HERE, f"{name}_{NAME}.csv")
        fields = []
        for rr in rows:
            for k in rr:
                if k not in fields:
                    fields.append(k)
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        print(f"[getResults] wrote {path} ({len(rows)} rows)")

    write("results", results)
    if history:
        write("history", history)
    if growth:
        write("growth", growth)

    groups = {}
    for r in results:
        groups.setdefault(r["row"], []).append(r)
    orders = []
    for row, recs in sorted(groups.items()):
        recs.sort(key=lambda r: r["N"])
        o = {"row": row, "N": " ".join(str(r["N"]) for r in recs)}
        for m in METRICS:
            vals = [r[m] for r in recs]
            o[f"{m}_values"] = " ".join(f"{v:.3e}" for v in vals)
            pw = []
            for k in range(len(recs) - 1):
                a, b = vals[k], vals[k + 1]
                hk, hk1 = recs[k]["h"], recs[k + 1]["h"]
                pw.append(math.log(a/b)/math.log(hk/hk1) if (a > 0 and b > 0) else float("nan"))
            o[f"p_{m}"] = " ".join(f"{p:.2f}" for p in pw)
        orders.append(o)
    write("orders", orders)

    print(f"[getResults] {'row':36s} N    L1(af)     p     L1(n)     p     L2(n)     p     L1(x)     p     L2(x)     p     E1        p     Emass     1-max     min       slv")
    for row, recs in sorted(groups.items()):
        for k, r in enumerate(recs):
            def pk(m):
                if k == 0:
                    return "  -- "
                a, b = recs[k - 1][m], r[m]
                return f"{math.log(a/b)/math.log(recs[k-1]['h']/r['h']):5.2f}" if (a > 0 and b > 0) else "  -- "
            print(f"[getResults] {row:36s} {r['N']:<4} {r['L1']:.3e} {pk('L1')} {r['LNormalDiff1']:.3e} {pk('LNormalDiff1')} "
                  f"{r['LNormalDiff2']:.3e} {pk('LNormalDiff2')} {r['LCentreDiff1']:.3e} {pk('LCentreDiff1')} "
                  f"{r['LCentreDiff2']:.3e} {pk('LCentreDiff2')} {r['E1']:.3e} {pk('E1')} {r['Emass']:.1e} "
                  f"{r.get('oneMinusMaxAlpha', float('nan')):.1e} {r.get('minAlpha', float('nan')):.1e} "
                  f"{str(r['nSliver_max'])}")
    if growth:
        print(f"[getResults] L1(alphaf) history (first {NHIST} steps) and the growth factors L1(k)/L1(k-1), davof rows:")
        for g in growth:
            l1 = [g[f"L1_{k}"] for k in range(1, NHIST + 1) if g[f"L1_{k}"] == g[f"L1_{k}"]]
            gs = [g[f"g_{k}"] for k in range(2, NHIST + 1) if g[f"g_{k}"] == g[f"g_{k}"]]
            print(f"[getResults] {g['row']} N{g['N']}: max growth {g['maxGrowth']:.2f}, "
                  f"geometric mean {g['geoMeanGrowth']:.2f}")
            print("[getResults]   L1: " + " ".join(f"{v:.1e}" for v in l1))
            print("[getResults]   g:  " + " ".join(f"{v:.2f}" for v in gs))
    return 0


if __name__ == "__main__":
    sys.exit(main())
