#!/usr/bin/env python3
"""Validate the exact piston reference (piston_reference.py) for this case.

1. On the axis, the rim integral equals the closed form.
2. Away from the baffle, it equals a direct two-dimensional quadrature of the
   Rayleigh integral over the piston surface.
3. Doubling the rim-integral nodes does not change the field in the physical
   region, including points next to the rim and the baffle.
4. The field is continuous across r = a, where the two forms meet.
"""

from __future__ import annotations

import math
import pathlib
import re
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

from piston_reference import piston_pressure, piston_pressure_on_axis

TOL = 1e-10


def case_params(path: pathlib.Path) -> dict[str, float]:
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        m = re.match(r"\s*([A-Z_]+)\s*=\s*([0-9eE+.\-]+)\s*$", line)
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def rayleigh_surface_quadrature(r, z, k, a, rho_c_u0, n_s=400, n_phi=800):
    xs, ws = leggauss(n_s)
    s, ws = 0.5 * a * (xs + 1.0), 0.5 * a * ws
    phi = np.linspace(0.0, 2.0 * math.pi, n_phi, endpoint=False)
    S, PH = np.meshgrid(s, phi, indexing="ij")
    out = np.empty(len(r), dtype=complex)
    for i, (ri, zi) in enumerate(zip(r, z)):
        R = np.sqrt(zi * zi + ri * ri + S * S - 2.0 * ri * S * np.cos(PH))
        out[i] = np.sum(np.exp(1j * k * R) / R * (S * ws[:, None])) * 2.0 * math.pi / n_phi
    # -(i omega rho / 2 pi) int e^{ikR}/R dS = -(i k / 2 pi) rho c u0 int ...
    return -1j * k * rho_c_u0 / (2.0 * math.pi) * out


def main() -> int:
    case = pathlib.Path(__file__).resolve().parent
    p = case_params(case / "caseParams.sh")
    k = 2.0 * math.pi * p["DRIVE_F"] / p["CG"]
    a = p["PISTON_RADIUS"]
    p0 = p["RHOG"] * p["CG"] * p["PISTON_U"]
    extent = p["PML_RMIN"]
    ok = True

    def check(name, err):
        nonlocal ok
        passed = err < TOL
        ok &= passed
        print(f"{'PASS' if passed else 'FAIL'}  {name}: max relative difference {err:.2e}")

    z = np.geomspace(1e-6, extent, 400)
    pe = piston_pressure(np.zeros_like(z), z, k, a, p0)
    check("rim integral vs on-axis closed form",
          np.max(np.abs(pe - piston_pressure_on_axis(z, k, a, p0))) / p0)

    rng = np.random.default_rng(1)
    r = np.r_[rng.uniform(0.0, extent, 30), 0.5 * a, 0.99 * a, 1.01 * a]
    zz = np.r_[rng.uniform(0.02, extent, 30), 0.05, 0.03, 0.03]
    pq = rayleigh_surface_quadrature(r, zz, k, a, p0)
    check("rim integral vs 2D Rayleigh quadrature (z >= 0.02 m)",
          np.max(np.abs(piston_pressure(r, zz, k, a, p0) - pq)) / p0)

    rg, zg = np.meshgrid(np.r_[0.0, np.linspace(1e-4, extent, 120), a - 1e-4, a + 1e-4],
                         np.r_[1e-4, np.linspace(1e-3, extent, 120)])
    p1 = piston_pressure(rg.ravel(), zg.ravel(), k, a, p0)
    p2 = piston_pressure(rg.ravel(), zg.ravel(), k, a, p0, n_phi=2048)
    check("node doubling 1024 -> 2048 (physical region, rim and baffle)",
          np.max(np.abs(p1 - p2)) / p0)

    # The gradient is singular at the rim, so the difference across a gap 2*eps*a
    # must vanish linearly with eps rather than be zero.
    zc = np.linspace(1e-3, extent, 50)
    jumps = []
    for eps in (1e-7, 1e-9, 1e-11):
        lo = piston_pressure(np.full_like(zc, a * (1 - eps)), zc, k, a, p0)
        hi = piston_pressure(np.full_like(zc, a * (1 + eps)), zc, k, a, p0)
        jumps.append(np.max(np.abs(lo - hi)) / p0)
    linear = all(j1 < 0.02 * j0 for j0, j1 in zip(jumps, jumps[1:]))
    ok &= linear
    print(f"{'PASS' if linear else 'FAIL'}  continuity across r = a: difference over the gap "
          f"a(1 -/+ eps) for eps = 1e-7, 1e-9, 1e-11: "
          + ", ".join(f"{j:.2e}" for j in jumps) + " (must fall with eps)")

    print("Reference validation", "PASSED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
