"""Exact pressure of a uniform circular piston in a rigid baffle.

Rayleigh integral, e^{-i omega t} convention:

    P(x) = -(i omega rho / 2 pi) int_piston u0 e^{ikR}/R dS.

In polar coordinates centred at the foot point of x on the baffle plane,
R dR = rho' drho', so the radial integral is exact and only a line integral
over the piston rim remains (plane wave minus edge wave).  With r the distance
of x from the axis and z its distance from the baffle:

  r <= a:  P = rho c u0 [ e^{ikz} - (1/2pi) int_0^{2pi} e^{ik R_e(phi)} dphi ],
           R_e = sqrt(z^2 + rho_e^2), rho_e = -r cos(phi) + sqrt(a^2 - r^2 sin^2(phi));

  r > a:   P = -(rho c u0 / 2pi) int_{-phi0}^{phi0}
               [ e^{ik R_2(phi)} - e^{ik R_1(phi)} ] dphi,   sin(phi0) = a/r,
           R_j = sqrt(z^2 + rho_j^2), rho_{1,2} = r cos(phi) -/+ sqrt(a^2 - r^2 sin^2(phi)).

On the axis, R_e = sqrt(z^2 + a^2) and P = rho c u0 [e^{ikz} - e^{ik sqrt(z^2+a^2)}].
The rim integrand is even in phi; it is integrated over [0, pi/2] and
[pi/2, pi] with Gauss-Legendre nodes, because rho_e develops a kink at
phi = pi/2 as r -> a.  The chord integral uses Gauss-Legendre nodes in t with
phi = phi0 sin(t), which removes the square-root endpoint behaviour.  Both
converge exponentially (checked by validate_reference.py).
"""

from __future__ import annotations

import functools
import math

import numpy as np
from scipy.special import roots_legendre


@functools.lru_cache(maxsize=8)
def _gauss_legendre(n: int) -> tuple[np.ndarray, np.ndarray]:
    return roots_legendre(n)


def piston_pressure(
    r: np.ndarray,
    z: np.ndarray,
    k: float,
    a: float,
    rho_c_u0: float,
    n_phi: int = 1024,
    chunk: int = 4096,
) -> np.ndarray:
    """Complex pressure amplitude at distance r from the axis and z from the baffle."""
    r = np.atleast_1d(np.asarray(r, dtype=float))
    z = np.atleast_1d(np.asarray(z, dtype=float))
    r, z = np.broadcast_arrays(r, z)
    out = np.empty(r.shape, dtype=complex)

    t, wt = _gauss_legendre(n_phi)
    # rim nodes on [0, pi/2] and [pi/2, pi]; weights include the 1/pi factor
    phi = np.r_[0.25 * math.pi * (t + 1.0), 0.25 * math.pi * (t + 3.0)]
    wphi = np.r_[wt, wt] / 4.0
    cphi, sphi = np.cos(phi), np.sin(phi)
    t, wt = 0.5 * math.pi * t, 0.5 * math.pi * wt
    st, ct = np.sin(t), np.cos(t)

    inside = r <= a
    idx_in = np.flatnonzero(inside)
    idx_out = np.flatnonzero(~inside)

    for s in range(0, idx_in.size, chunk):
        i = idx_in[s:s + chunk]
        ri, zi = r.flat[i][:, None], z.flat[i][:, None]
        rho_e = -ri * cphi + np.sqrt(np.maximum(a * a - (ri * sphi) ** 2, 0.0))
        edge = np.sum(wphi * np.exp(1j * k * np.sqrt(zi * zi + rho_e * rho_e)), axis=1)
        out.flat[i] = rho_c_u0 * (np.exp(1j * k * z.flat[i]) - edge)

    for s in range(0, idx_out.size, chunk):
        i = idx_out[s:s + chunk]
        ri, zi = r.flat[i][:, None], z.flat[i][:, None]
        phi0 = np.arcsin(a / ri)
        ph = phi0 * st
        d = np.sqrt(np.maximum(a * a - (ri * np.sin(ph)) ** 2, 0.0))
        rho1, rho2 = ri * np.cos(ph) - d, ri * np.cos(ph) + d
        f = np.exp(1j * k * np.sqrt(zi * zi + rho2 * rho2)) - np.exp(
            1j * k * np.sqrt(zi * zi + rho1 * rho1)
        )
        out.flat[i] = -rho_c_u0 / (2.0 * math.pi) * np.sum(wt * phi0 * ct * f, axis=1)

    return out


def piston_pressure_on_axis(z: np.ndarray, k: float, a: float, rho_c_u0: float) -> np.ndarray:
    """Closed form on the axis."""
    z = np.asarray(z, dtype=float)
    return rho_c_u0 * (np.exp(1j * k * z) - np.exp(1j * k * np.sqrt(z * z + a * a)))
