import gzip, os, re, sys
c = sys.argv[1]
def rd(p):
    if os.path.exists(p + ".gz"): return gzip.open(p + ".gz", "rt").read()
    return open(p).read()
def flist(txt, dtype=int):
    txt = txt[txt.index("FoamFile"):]
    txt = txt[txt.index("}") + 1:]
    m = re.search(r"\n\s*(\d+)\s*\n\s*\(", txt)
    n = int(m.group(1)); body = txt[m.end():]
    body = body[:body.index("\n)")]
    return n, [dtype(x) for x in body.split()]
hdr = rd(c + "/constant/polyMesh/owner")[:3000]
print("owner format:", re.search(r"format\s+(\w+)", hdr).group(1))
no, owner = flist(rd(c + "/constant/polyMesh/owner"))
nn, nei = flist(rd(c + "/constant/polyMesh/neighbour"))
cd = c + "/constant/cellDecomposition"
nc, proc = flist(rd(cd))
txt = rd(c + "/0/alpha.water")
i = txt.index("internalField")
m = re.search(r"List<scalar>\s*\n?\s*(\d+)\s*\n\s*\(", txt[i:])
n = int(m.group(1)); body = txt[i + m.end():]; body = body[:body.index("\n)")]
alpha = [float(x) for x in body.split()]
assert len(alpha) == nc == len(proc), (len(alpha), nc, len(proc))
iface = [1e-8 < a < 1 - 1e-8 for a in alpha]
onb = [False]*nc
for o, nb in zip(owner[:nn], nei):
    if proc[o] != proc[nb]: onb[o] = True; onb[nb] = True
nI = sum(iface); nIb = sum(1 for k in range(nc) if iface[k] and onb[k])
per = {}
for k in range(nc):
    if iface[k]: per[proc[k]] = per.get(proc[k], 0) + 1
print("cells %d, subdomains %d, cells per subdomain %s" % (nc, len(set(proc)), sorted(sum(1 for p in proc if p == q) for q in set(proc))))
print("interface cells (alpha in (1e-8, 1-1e-8)): %d, on %d subdomains %s" % (nI, len(per), dict(sorted(per.items()))))
print("interface cells with a face on a processor boundary: %d of %d; cells (any) on a processor boundary: %d" % (nIb, nI, sum(onb)))
