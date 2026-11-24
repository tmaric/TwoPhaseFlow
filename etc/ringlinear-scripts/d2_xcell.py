#!/usr/bin/env python3
"""the one cell of Part X (2026-11-24), read-only: python3 d2_xcell.py <serial E1> <parallel E1> <serial B> <parallel B> <time>"""
import re, sys
S1, P1, SB, PB_, TIME = sys.argv[1:6]
sys.argv = [sys.argv[0], S1, P1, TIME]
src = open("/work/scratch/tm83tomy/davof/agents/forkdefault/d2_ctl8-cmp.py").read()
exec(compile(src[:src.index("out = []\ndef w(")], "head", "exec"))
def uniform_value(txt, d):
    m_ = re.search(r"internalField\s+uniform\s+([^;]+);", txt)
    s_ = m_.group(1).replace("(", " ").replace(")", " ").split()
    return float(s_[0]) if d == 1 else np.array([float(x) for x in s_])
def pair(Sx, Px, name):
    global S, procs, cmap, cproc
    S = Sx
    plist = sorted(glob.glob(Px + "/processor*"), key=lambda s: int(s.rsplit("processor", 1)[1]))
    a, d, txta = internal(Sx + "/%s/%s" % (TIME, name))
    if a is None:
        u = uniform_value(txta, d); a = np.full(nC, u) if d == 1 else np.tile(u, (nC, 1))
    b = np.zeros_like(a)
    for ip, pd in enumerate(plist):
        v, dd, txtp = internal(pd + "/%s/%s" % (TIME, name))
        if v is None:
            u = uniform_value(txtp, dd); n_ = len(cmap[ip])
            v = np.full(n_, u) if dd == 1 else np.tile(u, (n_, 1))
        b[cmap[ip]] = v
    return a, b
h = 1.0/64
al = internal(S1 + "/%s/alpha.water" % TIME)[0]
al0 = internal(S1 + "/0/alpha.water" % () if False else S1 + "/0/alpha.water")[0]
m = internal(S1 + "/%s/areaNormalDavof.water" % TIME)[0]
for tag, (Sx, Px) in (("ringLinear 1 pass (E1)", (S1, P1)), ("ringSum (B)", (SB, PB_))):
    a, b = pair(Sx, Px, "interfaceCentre.water")
    d = np.sqrt(((a - b)**2).sum(1)); c = int(d.argmax())
    print("%s: maximum interfaceCentre difference %.3g (= %.3g h) in cell %d; alpha at time 0 %.17g, at the end %.17g; |m_c|/h^2 = %.6g; cell has a processor face: %s; rank %d; cells above 1e-12 h: %d" % (
        tag, d[c], d[c]/h, c, al0[c], al[c], np.linalg.norm(m[c])/h**2, "yes" if pb[c] else "no", cproc[c], (d > 1e-12*h).sum()))
    if tag.startswith("ringLinear"): cE1 = c
    else: cB = c
print("the cell of the ringLinear maximum is the cell of ringSum's maximum: %s (ringLinear cell %d, ringSum cell %d)" % ("yes" if cE1 == cB else "no", cE1, cB))
a, b = pair(SB, PB_, "interfaceCentre.water"); dB = np.sqrt(((a - b)**2).sum(1))
print("ringSum's difference in the ringLinear maximum cell: %.3g (= %.3g h)" % (dB[cE1], dB[cE1]/h))
def res(case):
    mx = []
    for l in open(case + "/log.advectorVoF", errors="replace"):
        g = re.search(r"residual max (\S+) mean", l)
        if g: mx.append(float(g.group(1)))
    return max(mx) if mx else None, len(mx)
print("positioningResidual: not written as a field; the log line of the positioning prints residual max: ringLinear 1 rank %s, 8 ranks %s (maximum over the run's lines, number of lines); ringSum 1 rank %s, 8 ranks %s" % (res(S1), res(P1), res(SB), res(PB_)))
