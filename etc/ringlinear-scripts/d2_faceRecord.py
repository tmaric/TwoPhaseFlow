#!/usr/bin/env python3
"""Part F (2026-11-24): the record of the processor faces whose end-time alphaf differs from serial.
python3 d2_faceRecord.py <serial case> <parallel case> <end time>. The helpers (parsing, mesh, decomposition) are those of d2_ctl8-cmp.py."""
import os, sys
src = open("/work/scratch/tm83tomy/davof/agents/forkdefault/d2_ctl8-cmp.py").read()
head = src[:src.index("out = []\ndef w(")]
exec(compile(head, "d2_ctl8-cmp-head", "exec"))
def sface(time, name):
    txt = rd(S + "/%s/%s" % (time, name)); i = txt.index("internalField")
    a = np.zeros(nF)
    if "nonuniform" in txt[i:i + 80]: a[:nInt] = readlist(txt, i, 1)
    else: a[:nInt] = float(re.search(r"internalField\s+uniform\s+([^;]+);", txt).group(1))
    bl = boundary_block(txt)
    for nm_, ty, nf, st in sb:
        v = patch_value(bl.get(nm_, ""), 1)
        if v is None: continue
        a[st:st + nf] = v[1] if v[0] == "list" else float(v[1])
    return a
def pfaces(time, name):
    vals = [[] for _ in range(nF)]
    for ip, pd in enumerate(procs):
        txt = rd(pd + "/%s/%s" % (time, name))
        pb_ = boundary_file(pd + "/constant/polyMesh/boundary")
        nIp = min(st for nm_, ty, nf, st in pb_); i = txt.index("internalField")
        if "nonuniform" in txt[i:i + 80]: iv = readlist(txt, i, 1)
        else: iv = np.full(nIp, float(re.search(r"internalField\s+uniform\s+([^;]+);", txt).group(1)))
        bl = boundary_block(txt); fm = fmap[ip]
        for k in range(nIp): vals[fm[k]].append((ip, iv[k]))
        for nm_, ty, nf, st in pb_:
            v = patch_value(bl.get(nm_, ""), 1)
            if v is None: continue
            arr = v[1] if v[0] == "list" else np.full(nf, float(v[1]))
            for k in range(nf): vals[fm[st + k]].append((ip, arr[k]))
    return vals
def cvec(root, procdirs, time, name):
    a, d, txt = internal(root + "/%s/%s" % (time, name))
    return a if a is not None else None
END = TIME
aS = sface(END, "alphaf.water"); aV = pfaces(END, "alphaf.water")
dif = [f for f in range(nInt) if aV[f] and max(abs(x[1] - aS[f]) for x in aV[f]) > 1e-12]
phiS = sface(END, "phi"); phiV = pfaces(END, "phi")
has_tag = os.path.exists(S + "/%s/alphafSource.water" % END) and os.path.exists(procs[0] + "/%s/alphafSource.water" % END)
if has_tag: tS = sface(END, "alphafSource.water"); tV = pfaces(END, "alphafSource.water")
al0 = internal(S + "/0/alpha.water")[0]; al1 = internal(S + "/%s/alpha.water" % END)[0]
nrm, _, _ = internal(S + "/%s/interfaceNormal.water" % END)
if nrm is None: nrm = np.zeros((nC, 3))
surf = np.sqrt((nrm**2).sum(1)) > 0
def sgn(x): return "0 (exactly zero)" if x == 0 else ("+" if x > 0 else "-")
print("%d processor faces differ from serial at the end time (above 1e-12); source tags on disk: %s" % (len(dif), "yes" if has_tag else "no"))
print()
for f in dif:
    o, n = owner[f], nei[f]
    print("face %d (serial internal face; owner cell %d on rank %d, neighbour cell %d on rank %d)" % (f, o, cproc[o], n, cproc[n]))
    print("- alphaf serial %.17g; sides: %s; the side equal to the serial value (1e-12): %s" % (aS[f], ", ".join("rank %d: %.17g" % (ip, x) for ip, x in aV[f]),
          ", ".join("rank %d" % ip for ip, x in aV[f] if abs(x - aS[f]) <= 1e-12) or "none"))
    print("- phi serial %.17g (%s); sides: %s" % (phiS[f], sgn(phiS[f]), ", ".join("rank %d: %.17g (%s)" % (ip, x, sgn(x)) for ip, x in phiV[f])))
    for tag, c in (("owner", o), ("neighbour", n)):
        print("- %s cell %d (rank %d): alpha before %.17g, after %.17g, surface cell of the step (written interfaceNormal nonzero): %s" % (tag, c, cproc[c], al0[c], al1[c], "yes" if surf[c] else "no"))
    if has_tag:
        print("- source tag: serial %g; sides: %s" % (tS[f], ", ".join("rank %d: %g" % (ip, x) for ip, x in tV[f])))
    print()
