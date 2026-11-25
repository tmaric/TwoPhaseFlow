#!/usr/bin/env python3
"""transport matrix with LIN = PURE's keys with ringLinear (2026-11-24): iso_lin/csv/iso_runs.csv + the runs -> ring_linear_transport_tables.md, PT1-PT4; error.dat columns 12, 13, 14 (corrected)"""
import csv, glob, math, re, statistics
W = "/work/scratch/tm83tomy/davof/agents/forkdefault"
OW = W + "/iso_lin"
def f(x):
    try: return float(x)
    except Exception: return float("nan")
M = {"RDF": "plicRDF", "PURE": "PURE", "LIN": "LIN (PURE keys, ringLinear, 1 pass)"}
FL = {"Shear2D": "shear2D", "Shear": "shear3D", "Deformation": "deformation3D"}
NS = {"Shear2D": (64, 128, 256), "Shear": (32, 64, 128), "Deformation": (32, 64, 128)}
RANKS = {("Shear2D", 64): 1, ("Shear2D", 128): 1, ("Shear2D", 256): 8, ("Shear", 32): 1, ("Shear", 64): 8, ("Shear", 128): 64,
         ("Deformation", 32): 1, ("Deformation", 64): 8, ("Deformation", 128): 64}
D = {}
for r in csv.DictReader(open(OW + "/csv/iso_runs.csv")):
    if r["method"] in M and r.get("mesh", "hex") == "hex":
        D.setdefault((r["flow"], int(r["N"]), int(r["co_index"]), r["method"]), []).append(r)
def first(k): return D[k][0] if k in D else None
def walls(k): return [f(r["wall_s"]) for r in D.get(k, [])]
def med(k): v = walls(k); return statistics.median(v) if v else float("nan")
def E(k): r = first(k); return f(r["Eshape_end"]) if r else float("nan")
def cdir(k): return f"{OW}/time/{k[0]}_hex_N{k[1]}_co{k[2]}/{k[3]}_r1"
def normals(k):
    try:
        e = [l.split() for l in open(glob.glob(cdir(k) + "/postProcessing/reconstructionError/*/error.dat")[0]) if not l.startswith("#")][-1]
        d1, d2, dinf = float(e[12]), float(e[13]), float(e[14])
        return d1, dinf, d2*d2/2, dinf*dinf/2
    except Exception:
        return None
DG = re.compile(r"^davof: ringLinear pass (\d+): y = 0 by solvability in (\d+) interface cells, negative weight in (\d+), max \|W \(q - qbar\)\.y\| (\S+)")
def diag(k):
    yz, ng, mw = [], [], []
    for l in open(cdir(k) + "/log.advectorVoF", errors="replace"):
        m = DG.match(l)
        if m: yz.append(int(m.group(2))); ng.append(int(m.group(3))); mw.append(float(m.group(4)))
    if not yz: return None
    return max(yz), yz[-1], max(ng), ng[-1], max(mw), mw[-1]
def order(a, b):
    if not (a == a and b == b) or a <= 0 or b <= 0: return ""
    return "%.2f" % math.log2(a/b)
def stamp():
    for c in glob.glob(f"{OW}/time/*/LIN_r1/.stamp"): return open(c).read().replace("\n", " ")
    return ""
out = []
def w(s=""): out.append(s)
w("## The transport matrix with LIN, the linear fit of the ring's area vectors (2026-11-24)")
w()
w("Cells: 3 flows x 3 N x OpenFOAM Co 0.2 and 0.5 (TARGET_ADVECTIVE_CO index 0 and 1), hexahedral meshes, the ranks of the curved matrix, one binary for every row and case (the build with `ringLinear`: %s), every davof key explicit, davofProcFaceRule symmetric, nAlphaBounds 6, clip false, snapTol 0. Rows: plicRDF; PURE; LIN = PURE's key set with normalModel ringLinear, ringLinearPasses 1, ringLinearCorrection true, ringLinearDiag true. One repetition per case, three for the six 3D N 128 cases on 64 ranks (median wall clock, all values listed). Through `d2_iso8.sbatch` (outputs `iso_lin/`). Normal norms from the corrected indices of error.dat (LNormalDiff1, LNormalDiff2, LNormalDiffInf = columns 13 to 15; 1 - n.n_e L1 = LNormalDiff2^2/2, Linf = LNormalDiffInf^2/2). Diagnostics of LIN (ringLinearDiag, per step): the interface cells with y = 0 by the solvability condition, the cells with a negative weight, the maximum of |W (q - qbar).y|; the maximum over the run and the value at the last step. Orders: log2 of the ratio between N and 2N." % stamp())
w()
w("| flow | N | ranks | Co label | solver max Courant | row | status | E_shape | order | L1(x) | order | abs(n-n_e) L1 | order | abs(n-n_e) Linf | 1-n.n_e L1 | 1-n.n_e Linf | mass error (%, max) | min alpha | max alpha - 1 | interface cells (last) | wall (s, median) | wall values (s) | wall ratio to plicRDF | E_shape ratio to plicRDF | E_shape ratio to PURE | interface cells LIN / PURE | y = 0 (max / last) | negative weight (max / last) | max abs(W (q-qbar).y) (max / last) |")
w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for fl in ("Shear2D", "Shear", "Deformation"):
    for ci in (0, 1):
        for m in ("RDF", "PURE", "LIN"):
            for i, N in enumerate(NS[fl]):
                k = (fl, N, ci, m); r = first(k)
                if r is None:
                    w("| %s | %d | %d | %s | | %s | not run |" % (FL[fl], N, RANKS[(fl, N)], "0.2" if ci == 0 else "0.5", M[m]) + " |"*22); continue
                nk = (fl, NS[fl][i+1], ci, m) if i + 1 < len(NS[fl]) else None
                nxt = first(nk) if nk else None
                nm = normals(k); nn = normals(nk) if nxt else None
                mw = med(k); wr = med((fl, N, ci, "RDF"))
                dg = diag(k) if m == "LIN" else None
                rp = first((fl, N, ci, "PURE"))
                w("| %s | %d | %d | %s | %.3f | %s | %s | %.4g | %s | %.4g | %s | %s | %s | %s | %s | %s | %.3g | %.3g | %.3g | %s | %.0f | %s | %.3f | %.3f | %.3f | %s | %s | %s | %s |" % (
                    FL[fl], N, RANKS[(fl, N)], r["co"], f(r["co_max_measured"]), M[m], r["status"], f(r["Eshape_end"]),
                    order(f(r["Eshape_end"]), f(nxt["Eshape_end"])) if nxt else "", f(r["L1x_end"]), order(f(r["L1x_end"]), f(nxt["L1x_end"])) if nxt else "",
                    ("%.3g" % nm[0]) if nm else "", order(nm[0], nn[0]) if (nm and nn) else "", ("%.3g" % nm[1]) if nm else "", ("%.3g" % nm[2]) if nm else "", ("%.3g" % nm[3]) if nm else "",
                    f(r["mass_dev_max_pct"]), f(r["min_alpha_run"]), f(r["max_excess_run"]), r["interface_cells_last"], mw, ", ".join("%.0f" % x for x in walls(k)), mw/wr,
                    E(k)/E((fl, N, ci, "RDF")), E(k)/E((fl, N, ci, "PURE")),
                    ("%.3f" % (f(r["interface_cells_last"])/f(rp["interface_cells_last"]))) if (m == "LIN" and rp) else "",
                    ("%d / %d" % (dg[0], dg[1])) if dg else "", ("%d / %d" % (dg[2], dg[3])) if dg else "", ("%.4g / %.4g" % (dg[4], dg[5])) if dg else ""))
w()
lab = lambda c: "%s N%d co%s" % (FL[c[0]], c[1], "0.2" if c[2] == 0 else "0.5")
cells = [(fl, N, ci) for fl in ("Shear2D", "Shear", "Deformation") for N in NS[fl] for ci in (0, 1)]
c3 = [c for c in cells if c[0] != "Shear2D"]
w("### The pre-registered statements")
w()
runs = [(c, first((*c, "LIN"))) for c in cells]
ok = []; bad = []
for c, r in runs:
    rp = first((*c, "PURE"))
    rat = f(r["interface_cells_last"])/f(rp["interface_cells_last"]) if r and rp else float("nan")
    good = r and r["status"] == "completed" and f(r["min_alpha_run"]) >= -1e-4 and f(r["max_excess_run"]) <= 1e-4 and abs(rat - 1) <= 0.10
    (ok if good else bad).append((c, r, rat))
w("- PT1 (every LIN run completes with min alpha >= -1e-4, max alpha - 1 <= 1e-4 and an interface cell count at the end within 10 %% of PURE's; falsifier: any run with more than 1.5 times PURE's interface cells, or min alpha below -1e-3): %d of %d runs satisfy all; exceptions: %s; runs with more than 1.5 times PURE's cells: %d; min alpha below -1e-3: %d." % (
  len(ok), len(cells), "; ".join("%s (%s, min %.3g, max-1 %.3g, cells %.3f x PURE)" % (lab(c), r["status"], f(r["min_alpha_run"]), f(r["max_excess_run"]), rat) for c, r, rat in bad) or "none",
  sum(1 for c, r, rat in ok + bad if rat > 1.5), sum(1 for c, r, rat in ok + bad if r and f(r["min_alpha_run"]) < -1e-3)))
lp = [(c, E((*c, "LIN"))/E((*c, "PURE"))) for c in c3]
sh = [(c, E((*c, "LIN"))/E((*c, "RDF"))) for c in c3 if c[0] == "Shear" and c[1] == 128]
w("- PT2 (LIN / PURE in E_shape below 1 in all twelve 3D cases; LIN / plicRDF on shear3D N 128 between 1.0 and 1.4 at both Co; falsifier of the second part: below 0.95 at both): LIN / PURE %s; below 1 in %d of 12. LIN / plicRDF at shear3D N 128: %s; between 1.0 and 1.4 at %d of 2; below 0.95 at %d of 2." % (
  "; ".join("%s %.3f" % (lab(c), x) for c, x in lp), sum(1 for c, x in lp if x < 1), "; ".join("%s %.3f" % (lab(c), x) for c, x in sh),
  sum(1 for c, x in sh if 1.0 <= x <= 1.4), sum(1 for c, x in sh if x < 0.95)))
big = [(c, med((*c, "LIN"))/med((*c, "RDF"))) for c in c3 if c[1] == 128]
w("- PT3 (LIN's median wall clock at N 128 on 64 ranks between 0.60 and 0.80 times plicRDF's in the four 3D cases): %s; between 0.60 and 0.80 in %d of 4." % ("; ".join("%s %.3f" % (lab(c), x) for c, x in big), sum(1 for c, x in big if 0.60 <= x <= 0.80)))
p4 = [(("Shear2D", N, ci), E(("Shear2D", N, ci, "LIN"))/E(("Shear2D", N, ci, "PURE"))) for N in (128, 256) for ci in (0, 1)]
w("- PT4 (shear2D LIN's E_shape at most PURE's at N 128 and 256, both Co): LIN / PURE %s; at most 1 in %d of 4." % ("; ".join("%s %.3f" % (lab(c), x) for c, x in p4), sum(1 for c, x in p4 if x <= 1)))
w()
w("### Courant numbers: the solver's own maximum against the target label")
w()
w("| flow | N | label 0.2 | label 0.5 |")
w("|---|---|---|---|")
for fl in ("Shear2D", "Shear", "Deformation"):
    for N in NS[fl]:
        a, b = first((fl, N, 0, "RDF")), first((fl, N, 1, "RDF"))
        w("| %s | %d | %.4f | %.4f |" % (FL[fl], N, f(a["co_max_measured"]) if a else float("nan"), f(b["co_max_measured"]) if b else float("nan")))
w()
open(W + "/ring_linear_transport_tables.md", "w").write("\n".join(out) + "\n")
print(len(out), "lines")
