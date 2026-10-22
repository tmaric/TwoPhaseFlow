#!/usr/bin/env python3
"""Where the unstable modes of an amplification matrix live.

Usage: localise_modes.py OUTDIR [N]

For the N eigenvalues of largest modulus outside the unit circle, prints
|g|, arg g, and the share of the scaled eigenvector in each state component
and in the physical domain and the PML cells with one, two or (in 3D)
three active rates.
"""
import json
import sys
from pathlib import Path

import numpy as np


def main():
    out = Path(sys.argv[1])
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    data = np.load(out / "G.npz")
    g, names, n_cells = data["G"], list(data["names"]), int(data["nCells"])
    s = json.loads((out / "setup.json").read_text())["setup"]
    h, c, n = s["h"], s["cG"], s["n"]

    scale = {"p": 1.0, "q": c / h, "psi": h / c, "pPrev": 1.0}
    d = np.concatenate([np.full(n_cells, scale.get(nm, 1.0 / h))
                        for nm in names])
    gs = g * (1 / d)[:, None] * d[None, :]
    values, vectors = np.linalg.eig(gs)
    order = np.argsort(-np.abs(values))
    order = [i for i in order if abs(values[i]) > 1 + 1e-10][:count]

    centres = -s["domainHalfWidth"] + (np.arange(n) + 0.5) * h
    in_pml = np.abs(centres) > s["halfWidth"]
    # blockMesh order: x fastest; count the active rates of each cell.
    grids = np.meshgrid(*[in_pml] * s["dim"], indexing="ij")
    rates = sum(g.astype(int) for g in grids).T.ravel()
    region = {"physical": rates == 0, "one rate": rates == 1,
              "two rates": rates == 2}
    if s["dim"] == 3:
        region["three rates"] = rates == 3

    for i in order:
        v = np.abs(vectors[:, i]) ** 2
        v /= v.sum()
        comp = {nm: v[k * n_cells:(k + 1) * n_cells].sum()
                for k, nm in enumerate(names)}
        per_cell = v.reshape(len(names), n_cells).sum(axis=0)
        where = {r: per_cell[m].sum() for r, m in region.items()}
        print(f"|g| {abs(values[i]):.6f}  arg {np.angle(values[i]):+.4f}  "
              + "  ".join(f"{k} {x:.2f}" for k, x in comp.items()) + "  |  "
              + "  ".join(f"{k} {x:.2f}" for k, x in where.items()))


if __name__ == "__main__":
    main()
