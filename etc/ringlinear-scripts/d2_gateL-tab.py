#!/usr/bin/env python3
"""ringLinear static gate (2026-11-24): gateL/<DIM>D_N<N>_p<POS>/{A,B,E1,E2} -> ring_linear_gate_tables.md, PL1-PL4, GATE; error.dat columns 12, 13, 14 (corrected)"""
import glob, math, os, re
W = "/work/scratch/tm83tomy/davof/agents/forkdefault"
G = W + "/gateL"; OLD = W + "/gate"
STATES = (("3D sphere", 3, (32, 64, 128)), ("2D circle", 2, (64, 128, 256, 512)))
ROWS = (("A", "plicRDF"), ("B", "ringSum (PURE)"), ("E1", "ringLinear, 1 pass"), ("E2", "ringLinear, 2 passes"), ("E0", "ringLinear, correction off (control)"))
out = []
def w(s=""): out.append(s)
def order(a, b):
    if a is None or b is None or a <= 0 or b <= 0: return ""
    return "%.2f" % math.log2(a/b)
NCD = re.compile(r"^NCD t (\S+) cell (\d+) err (\S+)")
DG = re.compile(r"^davof: ringLinear pass (\d+): y = 0 by solvability in (\d+) interface cells, negative weight in (\d+), max \|W \(q - qbar\)\.y\| (\S+)")
def lastline(c):
    fs = glob.glob(f"{c}/postProcessing/reconstructionError/*/error.dat")
    if not fs: return None
    ls = [l.split() for l in open(fs[0]) if l.strip() and not l.startswith("#")]
    return ls[-1] if ls else None
def read(c):
    r = lastline(c)
    if r is None or not os.path.exists(c + "/log.advectorVoF"): return None
    d = dict(d1=float(r[12]), d2=float(r[13]), dinf=float(r[14]), l1=float(r[2]), linf1=float(r[1]), nint=int(float(r[17])), raw=" ".join(r))
    d["c1"] = d["d2"]**2/2; d["cinf"] = d["dinf"]**2/2
    lg = open(c + "/log.advectorVoF", errors="replace").read()
    d["ok"] = "\nEnd" in lg
    ncd = {}
    for l in lg.splitlines():
        m = NCD.match(l)
        if m: ncd.setdefault(m.group(1), {})[int(m.group(2))] = float(m.group(3))
    d["n1e3"] = sum(1 for v in ncd[max(ncd, key=float)].values() if v > 1e-3) if ncd else None
    d["diag"] = [m.groups() for m in (DG.match(l) for l in lg.splitlines()) if m]
    return d
R = {}
for sn, dim, ns in STATES:
    for pos in (0, 1):
        for N in ns:
            for k, lab in ROWS:
                R[(dim, pos, N, k)] = read(f"{G}/{dim}D_N{N}_p{pos}/{k}")
w("## The area-weighted linear fit of the ring's area vectors, static gate (2026-11-24)")
w()
w("Route of the gate of 2026-11-23 (`d2_gateL.sbatch`): the exact initial state by initDavofFields, one step of deltaT/1000, reconstructionError, serial, one state per job, the rows one after the other, every davof key explicit, one binary for all rows (the build with `ringLinear` of commit 1638911: libVoF sha256 a4ef78aaa10a5b31, advectorVoF 616b35007350e317; also in etc/README-davof.md). Rows: plicRDF, ringSum (PURE), ringLinear with `ringLinearPasses 1` and 2 (`ringLinearCorrection true`, `ringLinearDiag true`). 3D sphere N 32/64/128; 2D circle N 64/128/256/512; standard and shifted centre ((0.013 0.007 0.003) in 3D; (0.013 0.007 0) in 2D, the y component lying along the circle's axis). Normal columns of error.dat with the corrected indices (LNormalDiff1, LNormalDiff2, LNormalDiffInf = columns 13 to 15; 1 - n.n_e L1 = LNormalDiff2^2/2, Linf = LNormalDiffInf^2/2). Cells above 1e-3: interface cells with 1 - n.n_e above 1e-3 (print-only key normalCellDiag). Orders: log2 of the ratio between N and 2N.")
w()
for sn, dim, ns in STATES:
    for pos in (0, 1):
        w(f"### {sn}, {'standard centre' if pos == 0 else 'shifted centre'}")
        w()
        w("| N | row | status | 1-n.n_e L1 | order | ratio to plicRDF | 1-n.n_e Linf | abs(n-n_e) L1 | order | ratio to plicRDF | abs(n-n_e) Linf | cells with 1-n.n_e > 1e-3 | interface cells | diagnostics (y = 0 by solvability / negative weight / max abs(W (q-qbar).y), per pass) |")
        w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for k, lab in ROWS:
            for i, N in enumerate(ns):
                d = R[(dim, pos, N, k)]
                if d is None and k == "E0":
                    continue
                if d is None:
                    w("| %d | %s | not run | | | | | | | | | | | |" % (N, lab)); continue
                a = R[(dim, pos, N, "A")]
                nx = R.get((dim, pos, ns[i+1], k)) if i + 1 < len(ns) else None
                o = lambda key: order(d[key], nx[key]) if nx else ""
                dg = "; ".join("pass %s: %s / %s / %s" % x for x in d["diag"][-(2 if k == "E2" else 1):]) if k.startswith("E") and d["diag"] else ""
                w("| %d | %s | %s | %.4g | %s | %.3f | %.4g | %.4g | %s | %.3f | %.4g | %s | %d | %s |" % (
                    N, lab, "completed" if d["ok"] else "no End", d["c1"], o("c1"), d["c1"]/a["c1"], d["cinf"], d["d1"], o("d1"), d["d1"]/a["d1"], d["dinf"],
                    d["n1e3"] if d["n1e3"] is not None else "", d["nint"], dg))
        w()
# the reproduction of the existing gate digits (A, B), and the checks
w("### Check: the rows plicRDF and ringSum against the gate table of 2026-11-23 (same states, serial)")
w()
same = 0; tot = 0; diff = []
for sn, dim, ns in STATES:
    for pos in (0, 1):
        for N in ns:
            for k in ("A", "B"):
                c = f"{OLD}/{dim}D_N{N}_p{pos}/{k}"
                o = lastline(c); d = R[(dim, pos, N, k)]
                if o is None or d is None: continue
                tot += 1
                if float(o[13]) == d["d2"] and float(o[12]) == d["d1"] and float(o[14]) == d["dinf"]: same += 1
                else: diff.append("%dD N%d p%d %s: %.6g against %.6g" % (dim, N, pos, k, float(o[13])**2/2, d["c1"]))
w("The last line of error.dat (LNormalDiff1, 2, Inf) of rows (a) and (b) of the new build is identical to the digit to the gate of 2026-11-23 in %d of %d states/rows; differences: %s." % (same, tot, "; ".join(diff) or "none"))
w()
w("### Check: 1 rank against 8 ranks, 3D sphere N 64, standard centre, ringLinear 1 pass (job 57155296)")
w()
def full(c):
    r = lastline(c); return [float(x) for x in r] if r else None
s1 = full(G + "/3D_N64_p0/E1"); s8 = full(G + "/par8_N64")
if s1 and s8:
    nm = {12: "LNormalDiff1 (abs(n-n_e) L1)", 13: "LNormalDiff2", 14: "LNormalDiffInf", 17: "nInterface"}
    w("| quantity | 1 rank | 8 ranks | relative difference |")
    w("|---|---|---|---|")
    mx = 0
    for i in (12, 13, 14, 17):
        rd = abs(s8[i] - s1[i])/abs(s1[i]) if s1[i] else float("nan"); mx = max(mx, rd)
        w("| %s | %.17g | %.17g | %.3g |" % (nm[i], s1[i], s8[i], rd))
    c1_, c8_ = s1[13]**2/2, s8[13]**2/2
    w("| 1-n.n_e L1 (LNormalDiff2^2/2) | %.17g | %.17g | %.3g |" % (c1_, c8_, abs(c8_ - c1_)/c1_))
    w()
    w("Agreement to 1e-12 relative: %s (largest relative difference %.3g). The columns LNormalInf_ and LNormal1_ of the 8-rank file are 0 (not reduced over the ranks)." % ("yes" if mx <= 1e-12 else "NO", mx))
else:
    w("not available")
w()
def ratio(dim, pos, N, k="E1", key="c1"):
    d, a = R.get((dim, pos, N, k)), R.get((dim, pos, N, "A"))
    return d[key]/a[key] if d and a else float("nan")
w("### The pre-registered statements")
w()
pl1 = [(pos, N, ratio(3, pos, N)) for pos in (0, 1) for N in (64, 128)]
w("- PL1 (ringLinear, 1 pass, 1-n.n_e L1 at most 0.5 times plicRDF's on the 3D sphere at N 64 and 128, both positions; falsifier: above 1.05 in any of the four): " + "; ".join("pos %d N %d: %.3f" % x for x in pl1) + ". At most 0.5 in %d of 4; above 1.05 in %d of 4." % (sum(1 for x in pl1 if x[2] <= 0.5), sum(1 for x in pl1 if x[2] > 1.05)))
pl2 = []
for pos in (0, 1):
    a, b = R[(3, pos, 64, "E1")], R[(3, pos, 128, "E1")]
    pl2.append((pos, math.log2(a["c1"]/b["c1"]) if a and b else float("nan")))
w("- PL2 (its order in 1-n.n_e L1 from N 64 to 128 at least 3, both positions; falsifier: below 2.5): " + "; ".join("pos %d: %.2f" % x for x in pl2) + ".")
pl3 = [(pos, N, ratio(2, pos, N)) for pos in (0, 1) for N in (128, 256)]
w("- PL3 (2D circle, ratio to plicRDF at most 0.5 at N 128 and 256, both positions; falsifier: above 1.05): " + "; ".join("pos %d N %d: %.3f" % x for x in pl3) + ". At most 0.5 in %d of 4; above 1.05 in %d of 4." % (sum(1 for x in pl3 if x[2] <= 0.5), sum(1 for x in pl3 if x[2] > 1.05)))
pl4 = []
for pos in (0, 1):
    for N in (32, 64, 128):
        a, b = R[(3, pos, N, "E1")], R[(3, pos, N, "E2")]
        if a and b: pl4.append((pos, N, b["c1"]/a["c1"] - 1))
w("- PL4 (two passes change the 3D L1 by less than 20 % against one pass): " + "; ".join("pos %d N %d: %+.1f %%" % (p, N, 100*x) for p, N, x in pl4) + ". Less than 20 %% in %d of %d." % (sum(1 for p, N, x in pl4 if abs(x) < 0.2), len(pl4)))
w()
w("Verdicts, literally against the written statements and falsifiers:")
w("- PL1: %s (at most 0.5 in %d of 4; falsifier, above 1.05 in any of the four: %s)." % ("held" if all(x[2] <= 0.5 for x in pl1) else "not held", sum(1 for x in pl1 if x[2] <= 0.5), "met" if any(x[2] > 1.05 for x in pl1) else "not met"))
w("- PL2: %s (order at least 3 at %d of 2 positions; falsifier, below 2.5: %s)." % ("held" if all(x[1] >= 3 for x in pl2) else "false", sum(1 for x in pl2 if x[1] >= 3), "met at %d of 2 positions" % sum(1 for x in pl2 if x[1] < 2.5) if any(x[1] < 2.5 for x in pl2) else "not met"))
w("- PL3: %s (at most 0.5 in %d of 4; falsifier, above 1.05: %s)." % ("held" if all(x[2] <= 0.5 for x in pl3) else "not held", sum(1 for x in pl3 if x[2] <= 0.5), "met" if any(x[2] > 1.05 for x in pl3) else "not met"))
w("- PL4: %s (less than 20 %% in %d of %d 3D states)." % ("held" if all(abs(x[2]) < 0.2 for x in pl4) else "false", sum(1 for x in pl4 if abs(x[2]) < 0.2), len(pl4)))
w()
w("### The GATE (at most 1.05 times plicRDF's 1-n.n_e L1 on the 3D sphere at N 64 and 128, both positions)")
w()
w("| candidate | pos 0 N 64 | pos 0 N 128 | pos 1 N 64 | pos 1 N 128 | passes |")
w("|---|---|---|---|---|---|")
for k, lab in ROWS[2:]:
    if k == "E0": continue
    rs = [ratio(3, pos, N, k) for pos in (0, 1) for N in (64, 128)]
    w("| %s | %s | %s |" % (lab, " | ".join("%.3f" % x for x in rs), "yes" if all(x <= 1.05 for x in rs) else "no"))
w()
w("GATE per point (ringLinear, 1 pass, ratio at most 1.05): " + "; ".join("pos %d N %d: %.3f %s" % (pos, N, ratio(3, pos, N, "E1"), "pass" if ratio(3, pos, N, "E1") <= 1.05 else "fail") for pos in (0, 1) for N in (64, 128)) + ".")
w()
open(W + "/ring_linear_gate_tables.md", "w").write("\n".join(out) + "\n")
print(len(out), "lines")
