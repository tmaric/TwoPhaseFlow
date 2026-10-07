"""Exact eigenmodes of a closed rigid gas/liquid cavity.

Gas on [0, xI], liquid on [xI, L]. With p_g = A cos(k_g x) and
p_l = B cos(k_l (L - x)), continuity of p and of (1/rho) dp/dx at xI gives

    Z_l sin(k_g L_g) cos(k_l L_l) + Z_g sin(k_l L_l) cos(k_g L_g) = 0,

the product form of Z_l tan(k_g L_g) + Z_g tan(k_l L_l) = 0.
"""
import math

import numpy as np


class TwoLayerMode:
    def __init__(self, rho_g, c_g, rho_l, c_l, x_i, length, omega, pa=1.0):
        self.rho_g, self.c_g, self.rho_l, self.c_l = rho_g, c_g, rho_l, c_l
        self.x_i, self.length, self.omega = x_i, length, omega
        self.k_g, self.k_l = omega / c_g, omega / c_l
        a = 1.0
        b = math.cos(self.k_g * x_i) / math.cos(self.k_l * (length - x_i))
        scale = pa / max(abs(a), abs(b))
        self.a, self.b = a * scale, b * scale

    def pressure(self, x):
        """Spatial shape; p = shape(x) cos(omega t)."""
        x = np.asarray(x, dtype=float)
        return np.where(
            x < self.x_i,
            self.a * np.cos(self.k_g * x),
            self.b * np.cos(self.k_l * (self.length - x)),
        )

    def velocity(self, x):
        """Spatial shape; U = shape(x) sin(omega t)."""
        x = np.asarray(x, dtype=float)
        return np.where(
            x < self.x_i,
            self.a * self.k_g / (self.rho_g * self.omega)
            * np.sin(self.k_g * x),
            -self.b * self.k_l / (self.rho_l * self.omega)
            * np.sin(self.k_l * (self.length - x)),
        )


def dispersion(omega, rho_g, c_g, rho_l, c_l, x_i, length):
    k_g, k_l = omega / c_g, omega / c_l
    l_g, l_l = x_i, length - x_i
    return (
        rho_l * c_l * math.sin(k_g * l_g) * math.cos(k_l * l_l)
        + rho_g * c_g * math.sin(k_l * l_l) * math.cos(k_g * l_g)
    )


def modes_near(f_target, rho_g, c_g, rho_l, c_l, x_i, length, band=0.2,
               samples=40001):
    """All eigen-angular-frequencies within f_target*(1 +- band)."""
    args = (rho_g, c_g, rho_l, c_l, x_i, length)
    w = np.linspace((1 - band), (1 + band), samples) * 2 * math.pi * f_target
    v = np.array([dispersion(wi, *args) for wi in w])
    roots = []
    for i in np.nonzero(v[:-1] * v[1:] < 0)[0]:
        lo, hi = w[i], w[i + 1]
        f_lo = v[i]
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            f_mid = dispersion(mid, *args)
            if f_lo * f_mid <= 0:
                hi = mid
            else:
                lo, f_lo = mid, f_mid
        roots.append(0.5 * (lo + hi))
    return roots


def select_mode(roots, f_target, which):
    w_target = 2 * math.pi * f_target
    below = [w for w in roots if w <= w_target]
    above = [w for w in roots if w > w_target]
    if which == "lower" and below:
        return max(below)
    if which == "upper" and above:
        return min(above)
    raise ValueError(f"no {which} mode near {f_target} Hz in {roots}")
