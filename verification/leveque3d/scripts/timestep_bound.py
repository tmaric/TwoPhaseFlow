#!/usr/bin/env python3
"""The exact maximum of the laboratory-frame speed, not a triangle bound.

    U = zdot + thetadot * (-r_z, 0, r_x) + Q(theta) * alpha(t) * u0(xi)

transcribed from movingFrameFlow.C, line for line. The triangle inequality
overstates this by about a third, because the rotation term peaks at the domain
corner where u0 is smallest, so the terms never reach their maxima together.

Sampling is legitimate here rather than a trial run: the field is given in
closed form and does not depend on the solution, so this is evaluation of a
known function, not measurement of a simulation.
"""

import math

import numpy as np

PI = np.pi


def u0_2d(x, z):
    inside = (x >= 0) & (x <= 1) & (z >= 0) & (z <= 1)
    sx, sz = np.sin(PI * x), np.sin(PI * z)
    return (np.where(inside, -sx * sx * np.sin(2 * PI * z), 0.0),
            np.where(inside, np.sin(2 * PI * x) * sz * sz, 0.0))


def yoff(t, amp, T):
    if amp == 0.0:
        return 0.0, 0.0
    s, w = t / T, 1.0 - t / T
    return amp * ((1 - s) * s + s * math.sin(w)), amp * ((1 - s) * s * s + s * math.sin(2 * w))


def yoffdot(t, amp, T):
    if amp == 0.0:
        return 0.0, 0.0
    s, w = t / T, 1.0 - t / T
    return ((amp / T) * (-s + (1 - s) + math.sin(w) - s * math.cos(w)),
            (amp / T) * (-s * s + (1 - s) * 2 * s + math.sin(2 * w) - 2 * s * math.cos(2 * w)))


def wmax_2d(T, period, revolutions, amp, centre=(0.5, 0.5), n=700, nt=800):
    g = (np.arange(n) + 0.5) / n
    X, Z = np.meshgrid(g, g, indexing="ij")
    thd = 2 * PI * revolutions / T
    best = 0.0
    for t in np.linspace(0.0, T, nt):
        th = 2 * PI * revolutions * t / T
        ct, st = math.cos(th), math.sin(th)
        a = math.cos(2 * PI * t / period) if period > 0 else 1.0
        yox, yoz = yoff(t, amp, T)
        ydx, ydz = yoffdot(t, amp, T)
        rx, rz = X - yox - centre[0], Z - yoz - centre[1]
        xix, xiz = centre[0] + (rx * ct + rz * st), centre[1] + (-rx * st + rz * ct)
        ubx, ubz = u0_2d(xix, xiz)
        ubx, ubz = a * ubx, a * ubz
        Ux = ydx + thd * (-rz) + (ubx * ct - ubz * st)
        Uz = ydz + thd * (rx) + (ubx * st + ubz * ct)
        best = max(best, float(np.max(np.hypot(Ux, Uz))))
    return best


def u0_3d(x, y, z, field, swirl=(0.5, 0.5)):
    if field == "leveque":
        Sx, Sy, Sz = (np.sin(2 * PI * v) for v in (x, y, z))
        sx, sy, sz = (np.sin(PI * v) ** 2 for v in (x, y, z))
        return 2 * sx * Sy * Sz, -Sx * sy * Sz, -Sx * Sy * sz
    inside = (x >= 0) & (x <= 1) & (y >= 0) & (y <= 1)
    sx, sy = np.sin(PI * x), np.sin(PI * y)
    ux = np.where(inside, np.sin(2 * PI * y) * sx * sx, 0.0)
    uy = np.where(inside, -np.sin(2 * PI * x) * sy * sy, 0.0)
    rho = np.sqrt((x - swirl[0]) ** 2 + (y - swirl[1]) ** 2)
    return ux, uy, (1.0 - 2.0 * rho) ** 2


def wmax_3d(T, period, revolutions, field, lz=1.0, centre=(0.5, 0.5), n=110, nt=240):
    g = (np.arange(n) + 0.5) / n
    X, Y, Z = np.meshgrid(g, g, g * lz, indexing="ij")
    thd = 2 * PI * revolutions / T
    best = 0.0
    for t in np.linspace(0.0, T, nt):
        th = 2 * PI * revolutions * t / T
        ct, st = math.cos(th), math.sin(th)
        a = math.cos(2 * PI * t / period) if period > 0 else 1.0
        rx, ry = X - centre[0], Y - centre[1]
        xix, xiy = centre[0] + (rx * ct + ry * st), centre[1] + (-rx * st + ry * ct)
        ux, uy, uz = u0_3d(xix, xiy, Z, field)
        ux, uy, uz = a * ux, a * uy, a * uz
        Ux = thd * (-ry) + (ux * ct - uy * st)
        Uy = thd * (rx) + (ux * st + uy * ct)
        best = max(best, float(np.max(np.sqrt(Ux ** 2 + Uy ** 2 + uz ** 2))))
    return best


CASES = [
    ("reversed-vortex", wmax_2d(8.0, 8.0, 1.0, 0.25), 8.0, (32, 64, 128, 256, 512, 1024, 2048)),
    ("deformation-sphere3d", wmax_3d(3.0, 6.0, 1.0, "spiralling"), 3.0, (32, 64, 128, 256)),
    ("leveque3d", wmax_3d(3.0, 6.0, 1.0, "leveque"), 3.0, (32, 64, 128, 256)),
]

print("exact max |w| over the domain and the period (framed configuration)\n")
for name, w, T, Ns in CASES:
    print(f"### {name}:  |w|max = {w:.6f}   (T = {T:g})")
    for N in Ns:
        dt_t = 0.2 / (N * w)
        steps = 4 * math.ceil(T / (4 * dt_t))
        dt = T / steps
        print(f"    N={N:5d}  steps={steps:8d}  dt={dt:.6e}  Co={w * dt * N:.4f}")
    print()
