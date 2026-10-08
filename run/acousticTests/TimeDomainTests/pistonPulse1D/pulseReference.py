"""Exact solution of a piston pulse through a planar interface.

Medium 1 (rho1, c1) fills 0 < x < xI, medium 2 (rho2, c2) the rest, which
is unbounded or ends in an ideal PML. The piston prescribes the velocity
U_p(t) and therefore reflects incident pressure like a rigid wall. With
R = (Z2 - Z1)/(Z2 + Z1), the pressure is a finite sum of travelling pulses:

  x < xI: Z1 sum_m R^m [U_p(t - (2m xI + x)/c1)
                        + R U_p(t - (2(m+1) xI - x)/c1)]
  x > xI: Z1 (1 + R) sum_m R^m U_p(t - (2m+1) xI/c1 - (x - xI)/c2)

The velocity is (right-going - left-going)/Z1 in medium 1 and p/Z2 in
medium 2; at x = 0 it equals U_p(t).
"""
import math

import numpy as np


class PistonPulse:
    def __init__(self, f0, n_cycles, u0):
        self.omega = 2.0 * math.pi * f0
        self.duration = n_cycles / f0
        self.u0 = u0

    def velocity(self, t):
        t = np.asarray(t, dtype=float)
        inside = (t > 0) & (t < self.duration)
        window = np.sin(math.pi * np.where(inside, t, 0) / self.duration) ** 4
        return np.where(inside, self.u0 * np.sin(self.omega * t) * window, 0.0)

    def acceleration_expression(self, density):
        """OpenFOAM expression for the piston boundary gradient.

        OpenFOAM prescribes dp/dn along the outward normal, here -x. With the
        outward normal velocity U_n = -U_p, dp/dn = -rho dU_n/dt = rho dU_p/dt.
        """
        w, tp, u0 = self.omega, self.duration, self.u0
        a = math.pi / tp
        return (
            f"(arg() < {tp!r}) ? {density!r}*{u0!r}*("
            f"{w!r}*cos({w!r}*arg())*pow(sin({a!r}*arg()),4)"
            f" + 4*{a!r}*sin({w!r}*arg())*pow(sin({a!r}*arg()),3)"
            f"*cos({a!r}*arg())) : 0"
        )


class LayeredPulseSolution:
    def __init__(self, pulse, rho1, c1, rho2, c2, x_interface):
        self.pulse = pulse
        self.rho1, self.c1, self.rho2, self.c2 = rho1, c1, rho2, c2
        self.z1, self.z2 = rho1 * c1, rho2 * c2
        self.x_i = x_interface
        self.r = (self.z2 - self.z1) / (self.z2 + self.z1)

    def _reflections(self, t_max):
        return int(math.ceil(t_max * self.c1 / (2.0 * self.x_i))) + 2

    def fields(self, x, t):
        """Pressure and velocity on the grid x (columns) at times t (rows)."""
        x = np.asarray(x, dtype=float)[None, :]
        t = np.asarray(t, dtype=float)[:, None]
        u = self.pulse.velocity
        r, xi, c1, c2 = self.r, self.x_i, self.c1, self.c2
        m_max = self._reflections(float(np.max(t)))
        right = np.zeros(np.broadcast(x, t).shape)
        left = np.zeros_like(right)
        trans = np.zeros_like(right)
        for m in range(m_max):
            rm = r**m
            right += rm * u(t - (2 * m * xi + x) / c1)
            left += rm * r * u(t - (2 * (m + 1) * xi - x) / c1)
            trans += rm * u(t - (2 * m + 1) * xi / c1 - (x - xi) / c2)
        medium1 = x < xi
        pressure = np.where(medium1, self.z1 * (right + left),
                            self.z1 * (1 + r) * trans)
        velocity = np.where(medium1, right - left,
                            self.z1 * (1 + r) * trans / self.z2)
        return pressure, velocity
