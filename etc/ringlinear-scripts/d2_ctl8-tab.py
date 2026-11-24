#!/usr/bin/env python3
"""controls of the 1 rank against 8 ranks difference of ringLinear (3D sphere N 64, standard centre): section appended to ring_linear_gate_tables.md"""
import glob, re, subprocess, sys
W = "/work/scratch/tm83tomy/davof/agents/forkdefault"
date = subprocess.check_output(["date", "+%F"]).decode().strip()
def last(c):
    fs = glob.glob(c + "/postProcessing/reconstructionError/*/error.dat")
    ls = [l.split() for l in open(fs[0]) if l.strip() and not l.startswith("#")]
    return [float(x) for x in ls[-1]]
def diag(c):
    out = None
    for l in open(c + "/log.advectorVoF", errors="replace"):
        m = re.match(r"davof: ringLinear pass (\d+): y = 0 by solvability in (\d+) interface cells, negative weight in (\d+), max \|W \(q - qbar\)\.y\| (\S+)", l)
        if m: out = m.groups()
    return out
G = W + "/gateL"
ROWS = (("H1: (a) plicRDF", G + "/3D_N64_p0/A", G + "/ctl8_A", False),
        ("H2: (b) ringSum", G + "/3D_N64_p0/B", G + "/ctl8_B", False),
        ("H3: ringLinear, correction off", W + "/gateL0/3D_N64_p0/E0", G + "/ctl8_E0", True),
        ("H4: ringLinear, 1 pass, correction on", G + "/3D_N64_p0/E1", G + "/par8_N64", True))
out = []
def w(s=""): out.append(s)
w()
w("## Controls of the 1-rank against 8-rank difference of ringLinear (%s)" % date)
w()
w("3D sphere N 64, standard centre, one step of deltaT/1000, binary libVoF sha256 a4ef78aaa10a5b31, advectorVoF 616b35007350e317. 8 ranks: the cases of the serial rows decomposed with `numberOfSubdomains 8; method scotch;` (the decomposeParDict of the case with no further coefficients; the same as job 57155296; scotch takes no n per direction), `mpirun -np 8 advectorVoF -parallel`; jobs 57168520 (A, B, E0 on 8 ranks), 57168521 (E0 on 1 rank), 57155296 (E1 on 8 ranks, reused). Only the reduced columns (LNormalDiff1, LNormalDiff2, LNormalDiffInf) and the interface-cell count are compared.")
w()
w("| row | quantity | 1 rank | 8 ranks | relative difference |")
w("|---|---|---|---|---|")
res = {}
names = {12: "LNormalDiff1", 13: "LNormalDiff2", 14: "LNormalDiffInf", 17: "nInterface"}
for lab, c1, c8, hasdiag in ROWS:
    a, b = last(c1), last(c8); mx = 0
    for i in (12, 13, 14, 17):
        rd = abs(b[i] - a[i])/abs(a[i]) if a[i] else float("nan"); mx = max(mx, rd)
        w("| %s | %s | %.17g | %.17g | %.3g |" % (lab, names[i], a[i], b[i], rd))
    res[lab] = mx
w()
w("Diagnostic lines of ringLinear (ringLinearDiag; last call of the run; cells with y = 0 by solvability / cells with a negative weight / max |W (q - qbar).y|):")
w()
w("| row | 1 rank | 8 ranks |")
w("|---|---|---|")
for lab, c1, c8, hasdiag in ROWS:
    if hasdiag:
        d1, d8 = diag(c1), diag(c8)
        f = lambda d: ("%s / %s / %s" % (d[1], d[2], d[3])) if d else "none"
        w("| %s | %s | %s |" % (lab, f(d1), f(d8)))
w()
w("Verdicts (agreement to 1e-12 relative in all four quantities):")
w()
for lab, c1, c8, hasdiag in ROWS:
    w("- %s: %s (largest relative difference %.3g)." % (lab, "agrees" if res[lab] <= 1e-12 else "does NOT agree", res[lab]))
w()
open(W + "/ctl8_section.md", "w").write("\n".join(out) + "\n")
print("\n".join(out))
