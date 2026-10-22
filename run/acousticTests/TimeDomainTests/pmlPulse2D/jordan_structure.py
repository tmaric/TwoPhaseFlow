#!/usr/bin/env python3
"""Jordan structure of an amplification matrix at g = 1.

Usage: jordan_structure.py OUTDIR [--chains]

On the active, scaled state of analyse_amplification.py, computes the
nullities of G - I and (G - I)^2 (relative tolerance 1e-12, with the gap
to the next singular value). A larger second nullity means Jordan chains:
G^n x = x + n (G - I) x for a generalized eigenvector x. With --chains,
also the size of (G - I) x for unit generalized vectors and where x and
(G - I) x live (components; PML cells by number of active rates).
Writes OUTDIR/jordan.json.
"""
import json
import sys
from pathlib import Path

import numpy as np

TOL = 1e-12


def active_scaled(out):
    data = np.load(out / "G.npz")
    g, names, n = data["G"], list(data["names"]), int(data["nCells"])
    s = json.loads((out / "setup.json").read_text())["setup"]
    h, c = s["h"], s["cG"]
    off = np.abs(g - np.diag(np.diag(g))).max(axis=1)
    active = ~((np.abs(np.diag(g) - 1) < 1e-14) & (off < 1e-14))
    scale = {"p": 1.0, "q": c / h, "psi": h / c, "pPrev": 1.0}
    d = np.concatenate([np.full(n, scale.get(nm, 1.0 / h)) for nm in names])
    gs = g * (1 / d)[:, None] * d[None, :]
    return gs[np.ix_(active, active)], active, names, n, s


def nullity(sv):
    rel = sv / sv.max()
    k = int(np.sum(rel < TOL))
    small = np.sort(rel)
    return k, float(small[k - 1]) if k else 0.0, float(small[k])


def shares(vectors, active, names, n, s):
    full = np.zeros((len(active), vectors.shape[1]))
    full[active] = vectors
    w = (full**2).sum(axis=1)
    w /= w.sum()
    comp = {nm: float(w[i * n:(i + 1) * n].sum()) for i, nm in enumerate(names)}
    cell = w.reshape(len(names), n).sum(axis=0)
    centres = -s["domainHalfWidth"] + (np.arange(s["n"]) + 0.5) * s["h"]
    in_pml = np.abs(centres) > s["halfWidth"]
    grids = np.meshgrid(*[in_pml] * s["dim"], indexing="ij")
    rates = sum(gr.astype(int) for gr in grids).T.ravel()
    return comp, {str(r): float(cell[rates == r].sum()) for r in range(s["dim"] + 1)}


def main():
    out = Path(sys.argv[1])
    chains = "--chains" in sys.argv
    gs, active, names, n, s = active_scaled(out)
    m = gs - np.eye(len(gs))
    result = {"size": int(active.sum())}
    if chains:
        _, s1, v1 = np.linalg.svd(m)
        _, s2, v2 = np.linalg.svd(m @ m)
    else:
        s1 = np.linalg.svd(m, compute_uv=False)
        s2 = np.linalg.svd(m @ m, compute_uv=False)
    k1, below1, above1 = nullity(s1)
    k2, below2, above2 = nullity(s2)
    result.update({"nullity1": k1, "gap1": [below1, above1],
                   "nullity2": k2, "gap2": [below2, above2],
                   "chains": k2 - k1})
    if chains and k2 > k1:
        n1, n2 = v1[-k1:].T, v2[-k2:].T
        u, _, _ = np.linalg.svd(n2 - n1 @ (n1.T @ n2), full_matrices=False)
        gen = u[:, :k2 - k1]
        image = m @ gen
        size = np.linalg.norm(image, axis=0)
        result["coupling"] = [float(size.min()), float(size.max())]
        result["generalized"] = dict(zip(("components", "rates"),
                                         shares(gen, active, names, n, s)))
        result["image"] = dict(zip(("components", "rates"),
                                   shares(image, active, names, n, s)))
    (out / "jordan.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"{out.name}: nullity G-I {k1}, (G-I)^2 {k2}, chains {k2 - k1}"
          + (f", coupling {result['coupling'][0]:.1e}..{result['coupling'][1]:.1e}"
             if "coupling" in result else ""))


if __name__ == "__main__":
    main()
