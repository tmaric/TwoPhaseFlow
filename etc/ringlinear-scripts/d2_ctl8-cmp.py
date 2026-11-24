#!/usr/bin/env python3
"""1 rank against 8 ranks, cell by cell and face by face (2026-11-24): python3 d2_ctl8-cmp.py <serial case> <parallel case> <time> -> markdown on stdout.
The processor fields are mapped to the serial mesh through cellProcAddressing and faceProcAddressing (read as ascii text, numpy); no reconstructPar."""
import glob, os, re, sys
import numpy as np
S, P, TIME = sys.argv[1], sys.argv[2], sys.argv[3]
T0 = "0"
TOL = 1e-12
def rd(p):
    return open(p).read()
def body_after(txt, start):
    m = re.compile(r"List<\w+>\s*\n?\s*(\d+)\s*\n?\s*\(").search(txt, start)
    n = int(m.group(1)); i = m.end()
    j = txt.index("\n)", i) if n > 0 else i
    return n, txt[i:j]
def readlist(txt, start, dim):
    n, b = body_after(txt, start)
    if n == 0: return np.zeros((0,) if dim == 1 else (0, dim))
    a = np.fromstring(b.replace("(", " ").replace(")", " "), sep=" ")
    return a if dim == 1 else a.reshape(-1, dim)
def dimof(txt):
    m = re.search(r"class\s+(\w+);", txt); c = m.group(1)
    return 3 if "Vector" in c else 1
def internal(path):
    txt = rd(path); d = dimof(txt)
    i = txt.index("internalField")
    seg = txt[i:i + 80]
    if "uniform" in seg and "nonuniform" not in seg:
        return None, d, txt
    return readlist(txt, i, d), d, txt
def boundary_block(txt):
    i = txt.index("boundaryField"); j = txt.index("{", i); depth = 0; k = j
    while True:
        if txt[k] == "{": depth += 1
        elif txt[k] == "}":
            depth -= 1
            if depth == 0: break
        k += 1
    blk = txt[j + 1:k]
    patches = {}
    pos = 0
    for m in re.finditer(r"\n\s*([\w\.]+)\s*\n\s*\{", blk):
        if m.start() < pos: continue
        s = m.end(); depth = 1; q = s
        while depth:
            if blk[q] == "{": depth += 1
            elif blk[q] == "}": depth -= 1
            q += 1
        patches[m.group(1)] = blk[s:q - 1]; pos = q
    return patches
def patch_value(body, d):
    m = re.search(r"\bvalue\s+(uniform|nonuniform)", body)
    if not m: return None
    if m.group(1) == "uniform":
        mm = re.search(r"value\s+uniform\s+([^;]+);", body)
        return ("uniform", mm.group(1))
    return ("list", readlist(body, m.start(), d))
def labels(path):
    txt = rd(path); i = txt.index("\n", txt.index("FoamFile"))
    m = re.search(r"\n\s*(\d+)\s*\n\s*\(", txt[i:])
    n = int(m.group(1)); b = txt[i + m.end(): i + txt[i + m.end():].index("\n)") + m.end()]
    return np.array(b.split(), dtype=np.int64)
def boundary_file(path):
    txt = rd(path); out = []
    for m in re.finditer(r"\n\s*([\w\.]+)\s*\n\s*\{([^}]*)\}", txt[txt.index("\n(", txt.index("FoamFile")):]):
        b = m.group(2)
        ty = re.search(r"\btype\s+(\w+);", b).group(1)
        nf = int(re.search(r"nFaces\s+(\d+);", b).group(1)); st = int(re.search(r"startFace\s+(\d+);", b).group(1))
        out.append((m.group(1), ty, nf, st))
    return out
# ---- the serial mesh
sm = S + "/constant/polyMesh/"
owner = labels(sm + "owner"); nei = labels(sm + "neighbour")
nInt = len(nei); nF = len(owner); nC = int(owner.max()) + 1
sb = boundary_file(sm + "boundary")
txt = rd(sm + "faces"); body = txt[txt.index("\n(", txt.index("FoamFile")) + 2:]
fp = np.array(re.findall(r"\d+\(([^)]*)\)", body), dtype=object)
faces = np.array([[int(x) for x in s.split()] for s in fp], dtype=np.int64)
assert faces.shape[1] == 4 and len(faces) == nF, faces.shape
# ---- the decomposition
procs = sorted(glob.glob(P + "/processor*"), key=lambda s: int(s.rsplit("processor", 1)[1]))
cproc = np.full(nC, -1, dtype=np.int64); cmap = []; fmap = []
for ip, pd in enumerate(procs):
    c = labels(pd + "/constant/polyMesh/cellProcAddressing"); cmap.append(c); cproc[c] = ip
    fmap.append(np.abs(labels(pd + "/constant/polyMesh/faceProcAddressing")) - 1)
assert (cproc >= 0).all()
isprocface = np.zeros(nF, dtype=bool)
isprocface[:nInt] = cproc[owner[:nInt]] != cproc[nei]
pb = np.zeros(nC, dtype=bool)
pb[owner[:nInt][isprocface[:nInt]]] = True; pb[nei[isprocface[:nInt]]] = True
# point-neighbour structure
cells_pairs = np.concatenate([np.stack([np.repeat(owner, 4), faces.ravel()], 1), np.stack([np.repeat(nei, 4), faces[:nInt].ravel()], 1)])
key = np.unique(cells_pairs[:, 0] * (faces.max() + 1) + cells_pairs[:, 1])
cell_of = key // (faces.max() + 1); pt_of = key % (faces.max() + 1)
ptmask = np.zeros(faces.max() + 1, dtype=bool); ptmask[pt_of[pb[cell_of]]] = True
ring_of_pb = np.zeros(nC, dtype=bool); ring_of_pb[cell_of[ptmask[pt_of]]] = True      # cells sharing a point with a processor-face cell (includes them)
ringonly = ring_of_pb & ~pb
out = []
def w(s=""): out.append(s)
w("Decomposition: %d subdomains; cells with a processor face: %d of %d; cells in the point-neighbour ring of those but without a processor face: %d; processor faces (serial internal faces between two subdomains): %d of %d internal faces." % (len(procs), pb.sum(), nC, ringonly.sum(), isprocface.sum(), nInt))
w()
def uniform_value(txt, d):
    m = re.search(r"internalField\s+uniform\s+([^;]+);", txt)
    s = m.group(1).replace("(", " ").replace(")", " ").split()
    return float(s[0]) if d == 1 else np.array([float(x) for x in s])
def cellfield(time, name):
    a, d, txt = internal(S + "/%s/%s" % (time, name))
    if a is None:
        u = uniform_value(txt, d); a = np.full(nC, u) if d == 1 else np.tile(u, (nC, 1))
    b = np.zeros_like(a)
    for ip, pd in enumerate(procs):
        v, dd, txtp = internal(pd + "/%s/%s" % (time, name))
        if v is None:
            u = uniform_value(txtp, dd); n = len(cmap[ip])
            v = np.full(n, u) if dd == 1 else np.tile(u, (n, 1))
        b[cmap[ip]] = v
    return a, b
def cellrow(label, time, name):
    ab = cellfield(time, name)
    if ab is None:
        w("| %s | %s | not written or uniform | | | | | |" % (label, time)); return None
    a, b = ab
    dif = np.abs(a - b) if a.ndim == 1 else np.sqrt(((a - b)**2).sum(1))
    m = dif > TOL
    w("| %s | %s | %d | %.3g | %d | %d | %d |" % (label, time, m.sum(), dif.max(), (m & pb).sum(), (m & ringonly).sum(), (m & ~ring_of_pb).sum()))
    return a, b
w("### Cell fields (serial against the processor fields mapped through cellProcAddressing; difference above 1e-12)")
w()
w("| field | time | cells differing | maximum difference | of those: with a processor face | only in the point-neighbour ring of one | elsewhere |")
w("|---|---|---|---|---|---|---|")
cellrow("alpha.water", T0, "alpha.water")
cellrow("alpha.water", TIME, "alpha.water")
res = {}
for nm in ("areaNormalDavof.water", "interfaceNormal.water", "interfaceCentre.water"):
    w("| %s | %s | not written at time 0 | | | | |" % (nm, T0))
    res[nm] = cellrow(nm, TIME, nm)
# ring sum recomputed offline from areaNormal (not the solver's field)
if res.get("areaNormalDavof.water"):
    for tag, k in (("serial", 0), ("8 ranks", 1)):
        pass
    alpha = internal(S + "/%s/alpha.water" % TIME)[0]
    cand = (alpha > 1e-8) & (alpha < 1 - 1e-8)
    def ring(m):
        d = np.zeros_like(m)
        # sum of m over the point-neighbour ring: cells sharing a point
        ptsum = np.zeros((faces.max() + 1, 3)); np.add.at(ptsum, pt_of, m[cell_of])
        # each point-neighbour counted once per shared point: use set semantics via the incidence of unique (cell, point) pairs
        return None
    ms, mp = res["areaNormalDavof.water"]
    # exact ring sum: for each candidate cell, the sum of m over the distinct cells sharing at least one point
    order = np.argsort(pt_of, kind="stable"); pt_sorted = pt_of[order]; cell_sorted = cell_of[order]
    starts = np.searchsorted(pt_sorted, np.arange(faces.max() + 2))
    cpts = {}
    idx = np.argsort(cell_of, kind="stable"); cs = cell_of[idx]; ps = pt_of[idx]; cst = np.searchsorted(cs, np.arange(nC + 1))
    cl = np.nonzero(cand)[0]
    ds = np.zeros((len(cl), 3)); dp = np.zeros((len(cl), 3))
    for r, c in enumerate(cl):
        nb = set()
        for p in ps[cst[c]:cst[c + 1]]:
            nb.update(cell_sorted[starts[p]:starts[p + 1]].tolist())
        nb = np.fromiter(nb, dtype=np.int64)
        ds[r] = ms[nb].sum(0); dp[r] = mp[nb].sum(0)
    dif = np.sqrt(((ds - dp)**2).sum(1)); m = dif > TOL
    cc = np.zeros(nC, dtype=bool); cc[cl[m]] = True
    w("| ringSum recomputed offline from the written areaNormal over the serial mesh's point-neighbour ring (candidate cells; NOT the solver's field) | %s | %d of %d | %.3g | %d | %d | %d |" % (TIME, m.sum(), len(cl), dif.max(), (cc & pb).sum(), (cc & ringonly).sum(), (cc & ~ring_of_pb).sum()))
w()
# ---- alphaf
w("### Face fractions alphaf.water (serial against the processor fields mapped through faceProcAddressing)")
w()
def serial_face(time):
    txt = rd(S + "/%s/alphaf.water" % time)
    i = txt.index("internalField")
    a = np.zeros(nF)
    if "nonuniform" in txt[i:i + 80]: a[:nInt] = readlist(txt, i, 1)
    else: a[:nInt] = float(re.search(r"internalField\s+uniform\s+([^;]+);", txt).group(1))
    bl = boundary_block(txt)
    for name, ty, nf, st in sb:
        v = patch_value(bl.get(name, ""), 1)
        if v is None: continue
        a[st:st + nf] = v[1] if v[0] == "list" else float(v[1])
    return a
def parallel_faces(time):
    # per face: the values from each side (processor faces appear in two subdomains); physical boundary and internal once
    vals = [[] for _ in range(nF)]
    for ip, pd in enumerate(procs):
        txt = rd(pd + "/%s/alphaf.water" % time)
        pb_ = boundary_file(pd + "/constant/polyMesh/boundary")
        nIp = min(st for nm, ty, nf, st in pb_)
        i = txt.index("internalField")
        if "nonuniform" in txt[i:i + 80]: iv = readlist(txt, i, 1)
        else: iv = np.full(nIp, float(re.search(r"internalField\s+uniform\s+([^;]+);", txt).group(1)))
        bl = boundary_block(txt); fm = fmap[ip]
        for k in range(nIp): vals[fm[k]].append((ip, iv[k]))
        for name, ty, nf, st in pb_:
            v = patch_value(bl.get(name, ""), 1)
            if v is None: continue
            arr = v[1] if v[0] == "list" else np.full(nf, float(v[1]))
            for k in range(nf): vals[fm[st + k]].append((ip, arr[k]))
    return vals
def alphaf_rows(time):
    a = serial_face(time); vals = parallel_faces(time)
    nsides = np.array([len(v) for v in vals])
    maxd = np.zeros(nF); sided = np.zeros(nF)
    for f in range(nF):
        v = vals[f]
        if not v: continue
        xs = [x[1] for x in v]
        maxd[f] = max(abs(x - a[f]) for x in xs)
        sided[f] = max(xs) - min(xs) if len(xs) > 1 else 0.0
    dif = maxd > TOL
    internal_rank = np.zeros(nF, dtype=bool); internal_rank[:nInt] = ~isprocface[:nInt]
    phys = np.zeros(nF, dtype=bool); phys[nInt:] = True
    w("| %s | %d | %d | %d | %d | %.3g | %d | %.3g |" % (time, dif.sum(), (dif & isprocface).sum(), (dif & internal_rank).sum(), (dif & phys).sum(),
      maxd[isprocface].max() if isprocface.any() else 0, (sided > TOL).sum(), sided.max()))
    return nsides, a, vals
w("| time | faces differing from serial (above 1e-12) | of those: processor faces | internal faces of a rank | physical boundary faces | max abs(delta alpha_f) on processor faces | processor faces whose two sides differ (above 1e-12) | max difference between the two sides |")
w("|---|---|---|---|---|---|---|---|")
ns0, _, _ = alphaf_rows(T0)
_, aE, valsE = alphaf_rows(TIME)
w()
w("Processor faces seen from two subdomains: %d of %d (each holds a value on both sides)." % ((ns0[isprocface] == 2).sum(), isprocface.sum()))

# ---- the area vectors of the state the reconstruction error is scored on (made after the advection of the step), recomputed offline
# from the written end-time face fractions: m_c = sum_f alpha_f S_f (S_f owner to neighbour, minus for the neighbour), each rank using its own side of a processor face
w()
w("### Area vectors recomputed offline from the end-time face fractions (the state scored by reconstructionError: the reconstruction made after the advection of the step)")
w()
txt = rd(sm + "points"); h = txt.index("}", txt.index("FoamFile")); mm_ = re.search(r"\n\s*(\d+)\s*\n\s*\(", txt[h:]); b = txt[h + mm_.end(): h + txt[h + mm_.end():].rindex("\n)") + mm_.end()]
pts = np.fromstring(b.replace("(", " ").replace(")", " "), sep=" ").reshape(-1, 3)
Sf = 0.5*np.cross(pts[faces[:, 2]] - pts[faces[:, 0]], pts[faces[:, 3]] - pts[faces[:, 1]])
def mvec(valfun):
    m = np.zeros((nC, 3))
    for f in range(nF):
        o = owner[f]
        m[o] += valfun(f, o)*Sf[f]
        if f < nInt:
            nn_ = nei[f]; m[nn_] -= valfun(f, nn_)*Sf[f]
    return m
ms = mvec(lambda f, c: aE[f])
def side_val(f, c):
    v = valsE[f]
    if len(v) == 1: return v[0][1]
    for ip, x in v:
        if ip == cproc[c]: return x
    return v[0][1]
mp = mvec(side_val)
dif = np.sqrt(((ms - mp)**2).sum(1)); m = dif > TOL
w("Area vectors m_c differing above 1e-12: %d cells (maximum difference %.3g); of those with a processor face: %d; only in the point-neighbour ring of one: %d; elsewhere: %d. Cells touching the differing face fractions (the faces' owner and neighbour cells): %d." % (
  m.sum(), dif.max(), (m & pb).sum(), (m & ringonly).sum(), (m & ~ring_of_pb).sum(),
  len(set(owner[[f for f in range(nF) if max(abs(x[1] - aE[f]) for x in valsE[f]) > TOL if valsE[f]]].tolist()) | set(nei[[f for f in range(nInt) if valsE[f] and max(abs(x[1] - aE[f]) for x in valsE[f]) > TOL]].tolist()))))
alpha = internal(S + "/%s/alpha.water" % TIME)[0]
cand = (alpha > 1e-8) & (alpha < 1 - 1e-8)
order = np.argsort(pt_of, kind="stable"); pt_sorted = pt_of[order]; cell_sorted = cell_of[order]
starts = np.searchsorted(pt_sorted, np.arange(faces.max() + 2))
idx = np.argsort(cell_of, kind="stable"); cs = cell_of[idx]; ps = pt_of[idx]; cst = np.searchsorted(cs, np.arange(nC + 1))
cl = np.nonzero(cand)[0]
ds = np.zeros((len(cl), 3)); dp = np.zeros((len(cl), 3))
for r, c in enumerate(cl):
    nb = set()
    for q in ps[cst[c]:cst[c + 1]]:
        nb.update(cell_sorted[starts[q]:starts[q + 1]].tolist())
    nb = np.fromiter(nb, dtype=np.int64)
    ds[r] = ms[nb].sum(0); dp[r] = mp[nb].sum(0)
dd = np.sqrt(((ds - dp)**2).sum(1)); mm = dd > TOL
cc = np.zeros(nC, dtype=bool); cc[cl[mm]] = True
w("Ring sums recomputed from these (candidate cells, alpha inside (1e-8, 1-1e-8)): differing above 1e-12: %d of %d cells (maximum difference %.3g); of those with a processor face: %d; only in the point-neighbour ring of one: %d; elsewhere: %d. Interface-cell normals (the ring sums' directions) differ in the same cells." % (
  mm.sum(), len(cl), dd.max(), (cc & pb).sum(), (cc & ringonly).sum(), (cc & ~ring_of_pb).sum()))
print("\n".join(out))
