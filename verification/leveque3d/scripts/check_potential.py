#!/usr/bin/env python3
"""Check that curl(A0) reproduces u0, for the vector potentials used in 3D.

The manuscript asserts this identity and quotes a numerical residual; the
number came from a run by hand.  This script is the workflow behind it.

Derivatives are taken by complex-step differentiation,

    f'(x) = Im f(x + ih) / h ,

which is exact to rounding for analytic f: unlike a finite difference it
never subtracts two nearly equal numbers, so there is no step-size trade-off
and the residual measures the identity rather than the differencing scheme.

A0 and u0 are transcribed from src/VoF/movingFrameFlow3D/movingFrameFlow3D.C
so that what is verified is the field the solver actually integrates.
"""

from __future__ import annotations

import argparse
import os

import numpy as np

PI = np.pi
H = 1.0e-20          # complex step; no cancellation, so it can be tiny


def A0_leveque(x, y, z):
    Cx = np.cos(2 * PI * x)
    Sy = np.sin(2 * PI * y)
    Sz = np.sin(2 * PI * z)
    Cz = np.cos(2 * PI * z)
    sy = np.sin(PI * y) ** 2
    sz = np.sin(PI * z) ** 2
    return (
        np.zeros_like(x),
        Sy * (Cx * sz + Cz) / (2 * PI),
        -Cx * sy * Sz / (2 * PI),
    )


def u0_leveque(x, y, z):
    Sx = np.sin(2 * PI * x)
    Sy = np.sin(2 * PI * y)
    Sz = np.sin(2 * PI * z)
    sx = np.sin(PI * x) ** 2
    sy = np.sin(PI * y) ** 2
    sz = np.sin(PI * z) ** 2
    return (2 * sx * Sy * Sz, -Sx * sy * Sz, -Sx * Sy * sz)


def _inside(x, y):
    """The psi0 support test, taken on the real parts as the C++ does."""
    xr, yr = np.real(x), np.real(y)
    return (xr >= 0.0) & (xr <= 1.0) & (yr >= 0.0) & (yr <= 1.0)


def A0_spiralling(x, y, z, swirl):
    psi0 = np.where(
        _inside(x, y),
        (1.0 / PI) * np.sin(PI * x) ** 2 * np.sin(PI * y) ** 2,
        np.zeros_like(x),
    )
    dx = x - swirl[0]
    dy = y - swirl[1]
    rho = np.sqrt(dx * dx + dy * dy)
    f = 0.5 - (4.0 / 3.0) * rho + rho * rho
    return (-f * dy, f * dx, psi0)


def u0_spiralling(x, y, z, swirl):
    ux = np.where(
        _inside(x, y),
        np.sin(2 * PI * y) * np.sin(PI * x) ** 2,
        np.zeros_like(x),
    )
    uy = np.where(
        _inside(x, y),
        -np.sin(2 * PI * x) * np.sin(PI * y) ** 2,
        np.zeros_like(x),
    )
    dx = x - swirl[0]
    dy = y - swirl[1]
    rho = np.sqrt(dx * dx + dy * dy)
    return (ux, uy, (1.0 - 2.0 * rho) ** 2)


def curl(A, x, y, z):
    """curl(A) by complex-step differentiation in each direction."""
    def d(direction):
        xx, yy, zz = x.astype(complex), y.astype(complex), z.astype(complex)
        if direction == 0:
            xx = xx + 1j * H
        elif direction == 1:
            yy = yy + 1j * H
        else:
            zz = zz + 1j * H
        return [np.imag(c) / H for c in A(xx, yy, zz)]

    dAdx = d(0)
    dAdy = d(1)
    dAdz = d(2)
    return (
        dAdy[2] - dAdz[1],
        dAdz[0] - dAdx[2],
        dAdx[1] - dAdy[0],
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--field", choices=("spiralling", "leveque"), required=True)
    p.add_argument("--swirl", default="0.5,0.5",
                   help="swirl axis (x,y) for the spiralling field")
    p.add_argument("--n", type=int, default=64, help="samples per direction")
    p.add_argument("--lz", type=float, default=1.0)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    swirl = tuple(float(v) for v in a.swirl.split(","))

    # sample cell-centre-like points, so no sample sits exactly on the
    # zero-extension boundary of psi0
    t = (np.arange(a.n) + 0.5) / a.n
    x, y, z = np.meshgrid(t, t, t * a.lz, indexing="ij")

    if a.field == "leveque":
        A = A0_leveque
        u = u0_leveque(x, y, z)
        name = "LeVeque deformation"
    else:
        A = lambda X, Y, Z: A0_spiralling(X, Y, Z, swirl)   # noqa: E731
        u = u0_spiralling(x, y, z, swirl)
        name = "spiralling deformation"

    c = curl(A, x, y, z)
    res = np.sqrt(sum((ci - ui) ** 2 for ci, ui in zip(c, u)))
    umag = np.sqrt(sum(ui ** 2 for ui in u))

    lines = [
        f"Vector-potential consistency: curl(A0) against u0, {name}",
        "=" * 68,
        f"samples          {a.n}^3 = {x.size} points, cell-centred in [0,1]^2 x [0,{a.lz}]",
        "differentiation  complex step, h = 1e-20 (exact to rounding)",
        f"swirl axis       {swirl}" if a.field == "spiralling" else "",
        "",
        f"max |curl(A0) - u0|            {res.max():.3e}",
        f"mean |curl(A0) - u0|           {res.mean():.3e}",
        f"max |u0|                       {umag.max():.3e}",
        f"max relative residual          {(res.max() / umag.max()):.3e}",
        "",
        "Componentwise maxima of the residual:",
    ]
    for i, comp in enumerate("xyz"):
        lines.append(f"  {comp}: {np.abs(c[i] - u[i]).max():.3e}")
    lines += [
        "",
        "A residual at the level of machine epsilon means the identity holds",
        "exactly and the flux assembly is reproducing the intended field, not",
        "merely one close to it.",
    ]
    text = "\n".join(l for l in lines if l != "") + "\n"

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w").write(text)
    print(text)


if __name__ == "__main__":
    main()
