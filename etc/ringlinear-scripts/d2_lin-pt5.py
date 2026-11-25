#!/usr/bin/env python3
"""PT5 (2026-11-24): shear3D N 64, Co 0.2 and 0.5, rows PURE and LIN, 8 ranks (iso_lin) against 1 rank (iso_lin1); error.dat columns 12, 13, 14"""
import csv, glob
W = "/work/scratch/tm83tomy/davof/agents/forkdefault"
def load(root):
    d = {}
    for r in csv.DictReader(open(root + "/csv/iso_runs.csv")):
        if r["flow"] == "Shear" and r["N"] == "64" and r["method"] in ("PURE", "LIN"):
            d[(int(r["co_index"]), r["method"])] = r
    return d
a8, a1 = load(W + "/iso_lin"), load(W + "/iso_lin1")
def l1(root, ci, m):
    e = [l.split() for l in open(glob.glob(f"{root}/time/Shear_hex_N64_co{ci}/{m}_r1/postProcessing/reconstructionError/*/error.dat")[0]) if not l.startswith("#")][-1]
    return float(e[13])**2/2
o = []
w = o.append
w("")
w("### PT5: shear3D N 64, 8 ranks against 1 rank (rows PURE and LIN; jobs 57171061, 57171062 for the 1-rank runs, 57171051, 57171052 for the 8-rank runs)")
w("")
w("| Co label | row | E_shape 8 ranks | E_shape 1 rank | ratio 8 / 1 | 1-n.n_e L1 8 ranks | 1-n.n_e L1 1 rank | ratio | min alpha 8 / 1 rank | wall 8 ranks (s) | wall 1 rank (s) |")
w("|---|---|---|---|---|---|---|---|---|---|---|")
bad = 0; rs = []
for ci in (0, 1):
    for m in ("PURE", "LIN"):
        r8, r1 = a8[(ci, m)], a1[(ci, m)]
        e8, e1 = float(r8["Eshape_end"]), float(r1["Eshape_end"])
        c8, c1 = l1(W + "/iso_lin", ci, m), l1(W + "/iso_lin1", ci, m)
        rs.append((ci, m, e8/e1)); bad += abs(e8/e1 - 1) > 0.01
        w("| %s | %s | %.6g | %.6g | %.4f | %.4g | %.4g | %.3f | %.3g / %.3g | %.0f | %.0f |" % ("0.2" if ci == 0 else "0.5", m, e8, e1, e8/e1, c8, c1, c8/c1, float(r8["min_alpha_run"]), float(r1["min_alpha_run"]), float(r8["wall_s"]), float(r1["wall_s"])))
w("")
w("PT5 (E_shape on 8 ranks within 1 %% of E_shape on 1 rank for both rows at both Co; falsifier: a difference above 1 %% in any of the four): ratios %s; differences above 1 %% in %d of 4: PT5 %s." % (
  "; ".join("Co %s %s %.4f" % ("0.2" if c == 0 else "0.5", m, x) for c, m, x in rs), bad, "held" if bad == 0 else "FALSE (falsifier met)"))
w("")
open(W + "/pt5_section.md", "w").write("\n".join(o) + "\n")
print("\n".join(o))
