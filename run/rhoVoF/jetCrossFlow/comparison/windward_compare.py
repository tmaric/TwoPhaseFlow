#!/usr/bin/env python3
"""The windward shapes of the simulated jets (windward_extract.py) against the
experiment (CXY_q7.51_calc.csv: X, Y, Z in metres, the windward shape of the
liquid jet in cross flow, ONERA; see README.md): one figure, and per run the
largest and the mean horizontal deviation |x_sim(y) - x_exp(y)| over the
experiment's heights that the simulated jet reaches (the paper: around 1 mm at
t = 7.1 ms on the fine mesh), the jet height of each curve.

    python3 windward_compare.py --exp CXY_q7.51_calc.csv RDF=windward_RDF.csv RING=windward_RING.csv --out windward.png [--title ...]
"""
import argparse, csv, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--exp", required=True)
ap.add_argument("runs", nargs="+", help="label=csv")
ap.add_argument("--out", default="windward.png")
ap.add_argument("--title", default="")
a = ap.parse_args()

def read_exp(p):
    pts = []
    for row in csv.reader(open(p)):
        if not row or row[0].startswith("#") or row[0].strip().lower() == "x": continue
        pts.append((float(row[0]), float(row[1])))
    pts.sort(key=lambda q: q[1])
    return np.array(pts)

def read_run(p):
    rows = list(csv.DictReader(open(p)))
    return np.array([(float(r["x_windward"]), float(r["y"])) for r in rows])

exp = read_exp(a.exp)
print("experiment: %d points, heights %.4f .. %.4f m" % (len(exp), exp[:, 1].min(), exp[:, 1].max()))
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(exp[:, 0]*1e3, exp[:, 1]*1e3, "k.", ms=4, label="experiment (windward shape)")
except ImportError:
    plt = None
for spec in a.runs:
    lab, p = spec.split("=", 1)
    run = read_run(p)
    if not len(run): print(lab, "empty"); continue
    ys = run[:, 1]; xs = run[:, 0]
    ok = (exp[:, 1] >= ys.min()) & (exp[:, 1] <= ys.max())
    xi = np.interp(exp[ok, 1], ys, xs)
    dev = np.abs(xi - exp[ok, 0])
    print("%-8s jet height %.4f m, %d of %d experiment heights covered, deviation max %.3e m mean %.3e m (x_sim - x_exp mean %+.3e m)" % (
        lab, ys.max(), ok.sum(), len(exp), dev.max() if len(dev) else float("nan"), dev.mean() if len(dev) else float("nan"), (xi - exp[ok, 0]).mean() if len(dev) else float("nan")))
    if plt: ax.plot(xs*1e3, ys*1e3, "-", lw=1.2, label=lab)
if plt:
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]"); ax.set_aspect("equal"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    if a.title: ax.set_title(a.title, fontsize=9)
    fig.tight_layout(); fig.savefig(a.out, dpi=150); print("wrote", a.out)
