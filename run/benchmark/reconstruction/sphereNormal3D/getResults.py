#!/usr/bin/env python3
"""Collect the sphere-normal reconstruction ladder into two CSVs and one figure.

Reads Cases/<scheme>/<setAlphaMethod>/N<N>/postProcessing/reconstructionError/0/error.dat
(the last row, columns by NAME from the `# Time ...` header line) and the case's
system/simulationParameter, and writes next to this script:

    results_sphereNormal3D.csv   scheme, setAlphaMethod, N, h, R_over_h, the 18 error.dat
                                 columns (LNormal1 ... nInterface), E_area_rel
    orders_sphereNormal3D.csv    least-squares + pairwise orders per (scheme, setAlphaMethod)
                                 for LNormalDiff1/2/Inf (|n - n_exact| mean/rms/max over the
                                 interface cells) and E_area_rel
    sphereNormal3D_L2.png        log-log LNormalDiff2 vs h

The normal error and the area are defined exactly as in leia's
leiaTestDavofNormal (E_L1_N / E_L2_N / E_LINF_N / E_AREA_REL), so
workflow/scripts/make_davof_normal_table.py --tpf-results results_sphereNormal3D.csv
merges these rows into the DAVOF order table.
"""
import csv
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def read_param(path):
    vals = {}
    for line in open(path):
        m = re.match(r"^\s*([A-Za-z_]\w*)\s+([^;]+);", line)
        if m:
            vals[m.group(1)] = m.group(2).strip()
    return vals


def read_error_dat(path):
    names, last = None, None
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            if s[1:].split() and s[1:].split()[0] == "Time":
                names = s[1:].split()
            continue
        last = s.split()
    if names is None or last is None:
        return None
    return dict(zip(names, last))


def fit(h, e):
    m = [(a, b) for a, b in zip(h, e) if a and b and a > 0 and b > 0]
    if len(m) < 2:
        return None, None
    lh, le = np.log([a for a, _ in m]), np.log([b for _, b in m])
    p = np.polyfit(lh, le, 1)
    r = le - np.polyval(p, lh)
    ss, st = float((r**2).sum()), float(((le - le.mean())**2).sum())
    return float(p[0]), (1.0 - ss/st if st > 0 else 1.0)


def pairwise(h, e):
    out = []
    for k in range(len(h) - 1):
        a, b = e[k], e[k + 1]
        out.append(np.log(a/b)/np.log(h[k]/h[k + 1]) if (a and b and a > 0 and b > 0) else None)
    return out


def main():
    rows = []
    for edat in sorted(glob.glob(os.path.join(HERE, "Cases", "*", "*", "*",
                                              "postProcessing", "reconstructionError", "0", "error.dat"))):
        case = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(edat))))
        par = read_param(os.path.join(case, "system", "simulationParameter"))
        rec = read_error_dat(edat)
        if rec is None:
            print(f"[getResults] no data in {edat}")
            continue
        n, L, R = int(par["N"]), float(par["L"]), float(par["R"])
        h = L/n
        row = {"scheme": par["RECONSCHEME"], "setAlphaMethod": par["SETALPHAMETHOD"],
               "N": n, "h": h, "R_over_h": R/h, "case": os.path.relpath(case, HERE)}
        for k, v in rec.items():
            row[k.rstrip("_")] = v
        areaSum, areaExact = float(rec.get("areaSum", 0)), float(rec.get("areaExact", 0))
        row["E_area_rel"] = abs(areaSum - areaExact)/areaExact if areaExact > 0 else ""
        rows.append(row)
    if not rows:
        print("[getResults] nothing found under Cases/"); return 1

    fields = list(rows[0])
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    rpath = os.path.join(HERE, "results_sphereNormal3D.csv")
    with open(rpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    print(f"[getResults] wrote {rpath} ({len(rows)} rows)")

    metrics = ["LNormalDiff1", "LNormalDiff2", "LNormalDiffInf", "E_area_rel"]
    groups = {}
    for r in rows:
        groups.setdefault((r["scheme"], r["setAlphaMethod"]), []).append(r)
    orders = []
    for (scheme, setA), recs in sorted(groups.items()):
        recs.sort(key=lambda r: -r["h"])
        h = [r["h"] for r in recs]
        o = {"scheme": scheme, "setAlphaMethod": setA, "N_coarsest": recs[0]["N"],
             "N_finest": recs[-1]["N"], "n_levels": len(recs),
             "nInterface_finest": recs[-1].get("nInterface", "")}
        for m in metrics:
            e = [float(r[m]) if r.get(m) not in ("", None) else None for r in recs]
            p, r2 = fit(h, e)
            o[f"{m}_finest"] = e[-1]
            o[f"p_{m}"] = p
            o[f"R2_{m}"] = r2
            o[f"pairwise_{m}"] = " ".join("--" if q is None else f"{q:.2f}" for q in pairwise(h, e))
        orders.append(o)
    opath = os.path.join(HERE, "orders_sphereNormal3D.csv")
    with open(opath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(orders[0]))
        w.writeheader(); w.writerows(orders)
    print(f"[getResults] wrote {opath}")
    for o in orders:
        print(f"[getResults] {o['scheme']:10s} {o['setAlphaMethod']:15s} N={o['N_coarsest']}..{o['N_finest']}"
              f"  p(L1)={o['p_LNormalDiff1']}  p(L2)={o['p_LNormalDiff2']}"
              f"  p(Linf)={o['p_LNormalDiffInf']}  p(A)={o['p_E_area_rel']}"
              f"  L2@finest={o['LNormalDiff2_finest']}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        colors = {"plicRDF": "#eb6834", "gradAlpha": "#1baf7a", "isoAlpha": "#eda100"}
        styles = {"cutCellIso": "-", "cutCellImpFunc": "--"}
        markers = {"plicRDF": "s", "gradAlpha": "D", "isoAlpha": "^"}
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        for (scheme, setA), recs in sorted(groups.items()):
            h = [r["h"] for r in recs]
            e = [float(r["LNormalDiff2"]) for r in recs]
            p, _ = fit(h, e)
            lab = f"{scheme}, {setA}" + (rf"  ($\propto h^{{{p:.2f}}}$)" if p is not None else "")
            ax.loglog(h, e, styles.get(setA, "-"), color=colors.get(scheme, "#e34948"),
                      marker=markers.get(scheme, "x"), ms=5, lw=1.6, label=lab)
        ax.set_xlabel(r"$h$ [m]"); ax.set_ylabel(r"$L_2$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$")
        ax.set_title("TwoPhaseFlow reconstruction of a sphere (R = 1 mm, 4R box)", fontsize=11)
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
        fig.tight_layout()
        fpath = os.path.join(HERE, "sphereNormal3D_L2.png")
        fig.savefig(fpath, dpi=200); plt.close(fig)
        print(f"[getResults] wrote {fpath}")
    except ImportError:
        print("[getResults] matplotlib not available; no figure")
    return 0


if __name__ == "__main__":
    sys.exit(main())
