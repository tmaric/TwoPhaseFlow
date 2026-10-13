#!/usr/bin/env python3
"""Spectral analysis of an amplification matrix formed by amplification.py.

Usage: analyse_amplification.py OUTDIR

Writes OUTDIR/spectrum.json with
- the spectral radius and the eigenvalues with |g| > 1 + 1e-10;
- the active state: components with an identity row (w where all PML
  rates vanish) stay at their initial value, zero for a physical initial
  state, and are removed; the remaining analysis uses the active state;
- the cluster at g = 1: algebraic multiplicity (|g - 1| < 1e-6) and the
  dimension of the null space of G - I; a difference means a Jordan block;
- max ||G^n|| over n <= 2^20 (n = 1..64 and powers of two), in the state
  scaled by p: 1, q: c/h, psi: h/c, w: 1/h, and the norms of the row
  blocks of the last power;
- for the theta rule, the semidiscrete eigenvalues
  lambda = (2/dt)(g - 1)/(g + 1) and their largest real part;
- for zero PML rates, the distance to the trapezoidal eigenvalues of the
  discrete Neumann Laplacian on the whole grid;
- the explicit-damping prediction max(1, |1 - sigma dt|) of the outermost
  PML cell.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np


def main():
    out = Path(sys.argv[1])
    data = np.load(out / "G.npz")
    g_mat, names = data["G"], list(data["names"])
    meta = json.loads((out / "setup.json").read_text())
    s = meta["setup"]
    n_cells, dt, h, c = meta["nCells"], s["deltaT"], s["h"], s["cG"]

    eig = np.linalg.eigvals(g_mat)
    mod = np.abs(eig)
    rho = float(mod.max())
    unstable = eig[mod > 1 + 1e-10]

    # Inert components: row j of G is the unit row e_j.
    size = len(g_mat)
    off = np.abs(g_mat - np.diag(np.diag(g_mat))).max(axis=1)
    inert = (np.abs(np.diag(g_mat) - 1) < 1e-14) & (off < 1e-14)
    active = ~inert
    g_act = g_mat[np.ix_(active, active)]
    eig_act = np.linalg.eigvals(g_act)

    near_one = int(np.sum(np.abs(eig_act - 1) < 1e-6))

    scale = {"p": 1.0, "q": c / h, "psi": h / c, "pPrev": 1.0}
    d = np.concatenate([np.full(n_cells, scale.get(nm, 1.0 / h))
                        for nm in names])[active]
    gs = g_act * (1 / d)[:, None] * d[None, :]
    # Null space of G - I in the scaled state, where all entries are O(1).
    sv = np.linalg.svd(gs - np.eye(len(gs)), compute_uv=False)
    null_dim = int(np.sum(sv < 1e-9 * sv.max()))
    norms = {}
    power = np.eye(len(gs))
    for n in range(1, 65):
        power = gs @ power
        norms[n] = float(np.linalg.norm(power, 2))
    n = 64
    # Stop squaring once the powers have grown beyond any bounded behaviour.
    while n < 2**20 and norms[n] < 1e100:
        power = power @ power
        n *= 2
        norms[n] = float(np.linalg.norm(power, 2))
    last = max(norms)
    worst = max(norms, key=norms.get)
    labels = np.concatenate([np.full(n_cells, i) for i in range(len(names))])[active]
    blocks = {names[i]: float(np.linalg.norm(power[labels == i], 2))
              for i in range(len(names)) if (labels == i).any()}

    result = {
        "method": s["method"], "theta": s["theta"], "courant": s["courant"],
        "sigmaMax": s["sigmaMax"], "deltaT": dt,
        "segregatedImplicitDamping": s.get("segregatedImplicitDamping", False),
        "size": size, "inert": int(inert.sum()),
        "superpositionError": meta["superpositionError"],
        "spectralRadius": rho,
        "nUnstable": int(len(unstable)),
        "unstable": [[float(z.real), float(z.imag)] for z in
                     sorted(unstable, key=abs, reverse=True)[:10]],
        "nearOne": {"algebraic": near_one, "geometric": null_dim},
        "powerNorm": {"max": norms[worst], "atN": worst,
                      "lastN": last, "last": norms[last], "blocksAtLast": blocks,
                      "history": {str(k): v for k, v in norms.items()}},
    }

    # Outermost PML cell: cubic profile at the cell centre.
    centres = -s["domainHalfWidth"] + (np.arange(s["n"]) + 0.5) * h
    depth = np.clip(np.abs(centres) - s["halfWidth"], 0, None)
    sigma_cell = (s["sigmaMax"] * (depth.max() / s["pmlThickness"]) ** s["pmlOrder"]
                  if s["pmlThickness"] > 0 else 0.0)
    result["outerCellSigmaDt"] = float(sigma_cell * dt)
    result["explicitDampingRadius"] = max(1.0, abs(1 - sigma_cell * dt))

    if s["method"] == "blockCoupled" and s["theta"] == 0.5:
        ok = np.abs(eig + 1) > 1e-12
        lam = 2 / dt * (eig[ok] - 1) / (eig[ok] + 1)
        result["semidiscreteMaxRe"] = float(lam.real.max())
        result["semidiscreteMaxReRelative"] = float(
            lam.real.max() / np.abs(lam).max())

    if s["sigmaMax"] == 0:
        m = np.arange(s["n"])
        k2 = (np.sin(np.pi * m / (2 * s["n"]))[:, None] ** 2
              + np.sin(np.pi * m / (2 * s["n"]))[None, :] ** 2)
        xi = (2 * c / h * np.sqrt(k2.ravel())) * dt
        xi = xi[xi > 0]
        ref = (1 + 0.5j * xi) / (1 - 0.5j * xi)
        ref = np.concatenate((ref, ref.conj()))
        dist = np.abs(ref[:, None] - eig[None, :]).min(axis=1)
        result["cavityEigenvalueError"] = float(dist.max())

    (out / "spectrum.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"{s['method']} Co {s['courant']:g} sigmaMax {s['sigmaMax']:g}"
          f"{' implicit damping' if result['segregatedImplicitDamping'] else ''}:"
          f" rho(G) - 1 = {rho - 1:.3e}, unstable {len(unstable)},"
          f" outer sigma dt {sigma_cell * dt:.3f}"
          f" (explicit-damping radius {result['explicitDampingRadius']:.4f})")
    print(f"  active {int(active.sum())} of {size}; g = 1: algebraic {near_one}, geometric {null_dim};"
          f" max ||G^n|| {norms[worst]:.3e} at n = {worst},"
          f" ||G^{last}|| {norms[last]:.3e}")
    print("  ||G^n||: " + ", ".join(f"{k}: {norms[k]:.2e}" for k in
                                     (1, 64, 1024, 16384, 2**20) if k in norms))
    print("  blocks of G^" + str(last) + ": "
          + ", ".join(f"{k} {v:.2e}" for k, v in blocks.items()))
    if "semidiscreteMaxRe" in result:
        print(f"  semidiscrete max Re lambda {result['semidiscreteMaxRe']:.3e}"
              f" (relative {result['semidiscreteMaxReRelative']:.2e})")
    if "cavityEigenvalueError" in result:
        print(f"  cavity eigenvalues: max distance {result['cavityEigenvalueError']:.2e}")


if __name__ == "__main__":
    main()
