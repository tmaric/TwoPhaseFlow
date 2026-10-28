#!/usr/bin/env python3
"""How far the material reaches: the bound the domain size rests on.

The manuscript sizes each domain by the furthest the material gets from the
rotation axis and, in 3D, how high it climbs. Those are properties of the exact
flow, so they are computed here from exact trajectories rather than measured
on a VOF field, where the alpha threshold and numerical wisps bias the extent.

The base flow is alpha(t) u0 with a scalar amplitude, so a body-frame
trajectory is the autonomous flow of u0 run to the pseudo-time
s(t) = int_0^t alpha = (P/2pi) sin(2pi t/P), which climbs to P/(2pi) at T/2 and
returns. Integrating u0 over [0, P/(2pi)] therefore visits every position the
material ever occupies. Distance to a z-parallel rotation axis is invariant
under the frame, so the body-frame extent is also the laboratory one.

The sphere is read from the case template's setAlphaFieldDict, so this bound
cannot drift from the configuration that is actually run. The maximum is taken
over the sphere's surface (the flow map is a homeomorphism) and is checked for
convergence in the number of surface seeds.
"""

from __future__ import annotations

import argparse
import math
import re

import numpy as np

PI = np.pi


def u0(p, field, swirl):
    x, y, z = p[:, 0], p[:, 1], p[:, 2]
    if field == "leveque":
        S = [np.sin(2 * PI * v) for v in (x, y, z)]
        s = [np.sin(PI * v) ** 2 for v in (x, y, z)]
        return np.stack([2 * s[0] * S[1] * S[2], -S[0] * s[1] * S[2],
                         -S[0] * S[1] * s[2]], 1)
    # 3D shear (Liovic et al. 2006): the 2D vortex horizontally, swirl axially
    inside = (x >= 0) & (x <= 1) & (y >= 0) & (y <= 1)
    ux = np.where(inside, np.sin(2 * PI * y) * np.sin(PI * x) ** 2, 0.0)
    uy = np.where(inside, -np.sin(2 * PI * x) * np.sin(PI * y) ** 2, 0.0)
    rho = np.hypot(x - swirl[0], y - swirl[1])
    return np.stack([ux, uy, (1 - 2 * rho) ** 2], 1)


def sphere(dict_path):
    s = open(dict_path).read()
    r = float(re.search(r"^\s*radius\s+([0-9.eE+-]+)\s*;", s, re.M).group(1))
    o = re.search(r"^\s*(?:origin|centre)\s+\(([^)]*)\)\s*;", s, re.M).group(1)
    return np.array([float(v) for v in o.split()]), r


def seeds(c, r, n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = PI * (1 + 5 ** 0.5) * i
    return np.stack([c[0] + r * np.cos(th) * np.sin(phi),
                     c[1] + r * np.sin(th) * np.sin(phi),
                     c[2] + r * np.cos(phi)], 1)


def extent(field, swirl, centre, radius, axis, period, n, nstep=600):
    p = seeds(centre, radius, n)
    h = period / (2 * PI) / nstep
    rmax = np.hypot(p[:, 0] - axis[0], p[:, 1] - axis[1]).max()
    lo, hi = p.min(0), p.max(0)
    for _ in range(nstep):
        k1 = u0(p, field, swirl)
        k2 = u0(p + 0.5 * h * k1, field, swirl)
        k3 = u0(p + 0.5 * h * k2, field, swirl)
        k4 = u0(p + h * k3, field, swirl)
        p = p + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        rmax = max(rmax, np.hypot(p[:, 0] - axis[0], p[:, 1] - axis[1]).max())
        lo, hi = np.minimum(lo, p.min(0)), np.maximum(hi, p.max(0))
    return float(rmax), lo, hi


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--field", required=True, choices=["leveque", "spiralling"])
    a.add_argument("--alpha-dict", required=True, help="template setAlphaFieldDict")
    a.add_argument("--axis", required=True, help='rotation centre, e.g. "(0.5 0.5 0)"')
    a.add_argument("--swirl", default="0.5 0.5")
    a.add_argument("--period", type=float, required=True)
    a.add_argument("--lz", type=float, required=True, help="domain height")
    a.add_argument("--resolutions", required=True, help="comma-separated N")
    a.add_argument("--out", required=True)
    args = a.parse_args()

    centre, radius = sphere(args.alpha_dict)
    ax = [float(v) for v in args.axis.strip("()").split()][:2]
    sw = [float(v) for v in args.swirl.replace(",", " ").split()][:2]
    Ns = [int(v) for v in args.resolutions.split(",")]

    rows = [(n, *extent(args.field, sw, centre, radius, ax, args.period, n))
            for n in (6000, 60000)]
    (n1, r1, lo1, hi1), (n2, r2, lo, hi) = rows
    # inscribed circle of the unit cross-section about the axis
    inscribed = min(ax[0], 1 - ax[0], ax[1], 1 - ax[1])
    margin = inscribed - r2
    top = args.lz - hi[2]

    with open(args.out, "w") as fh:
        w = fh.write
        w(f"# swept material extent, {args.field}; sphere centre {tuple(centre)} "
          f"radius {radius:g} (from {args.alpha_dict})\n")
        w(f"# rotation axis through ({ax[0]:g}, {ax[1]:g}), parallel to e_z; "
          f"period {args.period:g}\n")
        w(f"seeds {n1:6d}: r_max {r1:.4f}   z [{lo1[2]:.4f}, {hi1[2]:.4f}]\n")
        w(f"seeds {n2:6d}: r_max {r2:.4f}   z [{lo[2]:.4f}, {hi[2]:.4f}]\n")
        w(f"seed convergence: |dr_max| {abs(r2 - r1):.1e}, |dz_max| {abs(hi[2] - hi1[2]):.1e}\n")
        w(f"bounding box: x [{lo[0]:.4f}, {hi[0]:.4f}]  y [{lo[1]:.4f}, {hi[1]:.4f}]  "
          f"z [{lo[2]:.4f}, {hi[2]:.4f}]\n")
        w(f"horizontal margin to the inscribed circle (radius {inscribed:g}): {margin:.4f}\n")
        w(f"margin to the top of the domain (height {args.lz:g}): {top:.4f}\n")
        w("cells of margin:  N   horizontal   top\n")
        for N in Ns:
            w(f"               {N:4d}   {margin * N:8.1f}   {top * N:6.1f}\n")
    print(open(args.out).read())


if __name__ == "__main__":
    main()
