#!/usr/bin/env python3
"""Collect the reconstruction ladder of this benchmark into two CSVs and one figure.

Reads Cases/<scheme>/<setAlpha>/N<N>/postProcessing/reconstructionError/0/error.dat
(the last row, columns by NAME from the `# Time ...` header line), the case's
system/simulationParameter and, for the davof scheme, the explicit-plane
defects printed by the scheme in log.reconstructInterface, and writes next to
this script (NAME = the directory name, sphereNormal3D or planeNormal3D):

    results_<NAME>.csv   scheme, setAlpha, N, h, R_over_h, the error.dat columns
                         (LNormal1 ... LCentreDiffInf), E_area_rel,
                         maxPlaneVolDefect, maxPlaneAreaDefect (davof rows)
    orders_<NAME>.csv    least-squares + pairwise orders per (scheme, setAlpha)
                         for LNormalDiff1/2/Inf (|n - n_exact| mean/rms/max over
                         the interface cells), LCentreDiff1/2/Inf (the distance of
                         the reconstructed centre from the surface) and E_area_rel
    <NAME>_L2.png        log-log LNormalDiff2 and LCentreDiff2 vs h

planeNormal3D is the round-off gate: every davof row must have LNormalDiffInf,
LCentreDiff2/h and the two explicit-plane defects below GATE (1e-12) and
LCentreDiffInf/h below GATE_MAX (1e-10), or the script exits non-zero. The
looser maximum is the floor of the alpha-positioned plane in corner-sliver
cells (measured 2026-10-12 at N = 20: three cells with 1 - alpha = 1.8e-6 and a
polygon of 5e-4 h^2 set 1.3e-11 h; the L2 is 7.7e-13 h, the median 1e-13 h;
leia records the same 1e-11 h floor for its explicit position).

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
NAME = os.path.basename(HERE)
GATE = 1e-12
GATE_MAX = 1e-10
METRICS = ["LNormalDiff1", "LNormalDiff2", "LNormalDiffInf",
           "LCentreDiff1", "LCentreDiff2", "LCentreDiffInf", "E_area_rel"]


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


def read_davof_defects(log):
    """The last 'davof: ... explicit plane defects ...' line of the log."""
    vol = area = None
    if not os.path.isfile(log):
        return vol, area
    for line in open(log):
        m = re.search(r"max \|alpha_plane - alpha\| = ([-+.eE0-9]+), max \|A n - m\|/V\^\(2/3\) = ([-+.eE0-9]+)", line)
        if m:
            vol, area = float(m.group(1)), float(m.group(2))
    return vol, area


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


def set_alpha_label(par):
    if par["SETALPHAMETHOD"] == "davof":
        return "davof-" + par["DAVOFMETHOD"]
    return par["SETALPHAMETHOD"]


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
        n, L = int(par["N"]), float(par["L"])
        R = float(par["R"]) if "R" in par else None
        h = L/n
        row = {"scheme": par["RECONSCHEME"], "setAlphaMethod": set_alpha_label(par),
               "N": n, "h": h, "R_over_h": (R/h if R else ""), "case": os.path.relpath(case, HERE)}
        for k, v in rec.items():
            row[k.rstrip("_")] = v
        areaSum, areaExact = float(rec.get("areaSum", 0)), float(rec.get("areaExact", 0))
        row["E_area_rel"] = abs(areaSum - areaExact)/areaExact if areaExact > 0 else ""
        vol, area = (None, None)
        if par["RECONSCHEME"] == "davof":
            vol, area = read_davof_defects(os.path.join(case, "log.reconstructInterface"))
        row["maxPlaneVolDefect"] = "" if vol is None else vol
        row["maxPlaneAreaDefect"] = "" if area is None else area
        rows.append(row)
    if not rows:
        print("[getResults] nothing found under Cases/"); return 1

    fields = list(rows[0])
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    rpath = os.path.join(HERE, f"results_{NAME}.csv")
    with open(rpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    print(f"[getResults] wrote {rpath} ({len(rows)} rows)")

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
        for m in METRICS:
            e = [float(r[m]) if r.get(m) not in ("", None) else None for r in recs]
            p, r2 = fit(h, e)
            o[f"{m}_finest"] = e[-1]
            o[f"p_{m}"] = p
            o[f"R2_{m}"] = r2
            o[f"pairwise_{m}"] = " ".join("--" if q is None else f"{q:.2f}" for q in pairwise(h, e))
        orders.append(o)
    opath = os.path.join(HERE, f"orders_{NAME}.csv")
    with open(opath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(orders[0]))
        w.writeheader(); w.writerows(orders)
    print(f"[getResults] wrote {opath}")

    def fmt(p):
        return "  -- " if p is None else f"{p:5.2f}"
    print(f"[getResults] {'scheme':10s} {'setAlpha':22s} N        p(n L1) p(n L2) p(n Linf) p(x L1) p(x L2) p(x Linf)"
          f"  L2(n)@finest  L2(x)@finest  Linf(x)/h@finest")
    for o in orders:
        hf = float(o["N_finest"])
        hfin = [r["h"] for r in groups[(o["scheme"], o["setAlphaMethod"])] if r["N"] == o["N_finest"]][0]
        print(f"[getResults] {o['scheme']:10s} {o['setAlphaMethod']:22s} {o['N_coarsest']:>3}..{o['N_finest']:<4}"
              f"  {fmt(o['p_LNormalDiff1'])}  {fmt(o['p_LNormalDiff2'])}  {fmt(o['p_LNormalDiffInf'])}"
              f"    {fmt(o['p_LCentreDiff1'])}  {fmt(o['p_LCentreDiff2'])}  {fmt(o['p_LCentreDiffInf'])}"
              f"    {o['LNormalDiff2_finest']:.3e}     {o['LCentreDiff2_finest']:.3e}     "
              f"{o['LCentreDiffInf_finest']/hfin:.3e}")

    status = 0
    if NAME == "planeNormal3D":
        print(f"[getResults] round-off gate (GATE = {GATE:g}, GATE_MAX = {GATE_MAX:g} on Linf(dist)/h): davof rows")
        for r in rows:
            if r["scheme"] != "davof":
                continue
            nerr = float(r["LNormalDiffInf"])
            x2 = float(r["LCentreDiff2"])/r["h"]
            xerr = float(r["LCentreDiffInf"])/r["h"]
            vol = float(r["maxPlaneVolDefect"]) if r["maxPlaneVolDefect"] != "" else float("nan")
            area = float(r["maxPlaneAreaDefect"]) if r["maxPlaneAreaDefect"] != "" else float("nan")
            ok = nerr < GATE and x2 < GATE and xerr < GATE_MAX and vol < GATE and area < GATE
            status |= (0 if ok else 2)
            print(f"[getResults]   {r['setAlphaMethod']:22s} N={r['N']:<4} Linf|n - n_exact| = {nerr:.2e}"
                  f"  L2(dist)/h = {x2:.2e}  Linf(dist)/h = {xerr:.2e}  max|alpha_plane - alpha| = {vol:.2e}"
                  f"  max|A n - m|/V^(2/3) = {area:.2e}  {'PASS' if ok else 'FAIL'}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        colors = {"plicRDF": "#eb6834", "gradAlpha": "#1baf7a", "isoAlpha": "#eda100", "davof": "#2a6ebb"}
        styles = {"cutCellIso": "-", "cutCellImpFunc": "--", "davof-detrixheAslam": "-",
                  "davof-quadraticFaces": "-.", "davof-exactSphere": ":"}
        markers = {"plicRDF": "s", "gradAlpha": "D", "isoAlpha": "^", "davof": "o"}
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
        for ax, metric, ylab in ((axes[0], "LNormalDiff2", r"$L_2$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$"),
                                 (axes[1], "LCentreDiff2", r"$L_2$ of the centre distance [m]")):
            for (scheme, setA), recs in sorted(groups.items()):
                h = [r["h"] for r in recs]
                e = [float(r[metric]) for r in recs]
                p, _ = fit(h, e)
                lab = f"{scheme}, {setA}" + (rf"  ($\propto h^{{{p:.2f}}}$)" if p is not None else "")
                ax.loglog(h, e, styles.get(setA, "-"), color=colors.get(scheme, "#e34948"),
                          marker=markers.get(scheme, "x"), ms=5, lw=1.6, label=lab)
            ax.set_xlabel(r"$h$ [m]"); ax.set_ylabel(ylab)
            ax.grid(True, which="both", ls=":", alpha=0.4)
        axes[0].set_title(f"TwoPhaseFlow reconstruction: {NAME}", fontsize=11)
        axes[1].legend(frameon=False, fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5))
        fig.tight_layout()
        fpath = os.path.join(HERE, f"{NAME}_L2.png")
        fig.savefig(fpath, dpi=200); plt.close(fig)
        print(f"[getResults] wrote {fpath}")
    except ImportError:
        print("[getResults] matplotlib not available; no figure")
    return status


if __name__ == "__main__":
    sys.exit(main())
