#!/usr/bin/env python3
"""max ||w||_1 of the laboratory field: the quantity that sets the Courant number.

OpenFOAM forms the interface Courant number as

    Co = 0.5 * sum_f |phi_f| * dt / V ,

and on a Cartesian cell each direction contributes two faces of area A with
V = A*h, so Co = (|w_x| + |w_y| + |w_z|) * dt / h. It is the 1-norm of the
velocity that matters, not its magnitude -- sizing the step from the magnitude
leaves the reported Courant number a factor ||w||_1/||w||_2 too high, which is
1.31 for the 2D case and 1.47 for the 3D shear.

Taking the maximum over the whole domain rather than over the interface keeps
the step independent of the solution: it is a property of the case definition,
computable before anything is run, and it bounds the Courant number in every
cell rather than only where the interface happens to be.

    dt = Co * h / max||w||_1 .

u0 and the frame kinematics are transcribed from movingFrameFlow.C and
movingFrameFlow3D.C, so this bounds the field the solver integrates.
"""

import math

import numpy as np

PI = np.pi


def u0_2d(x, z):
    inside = (x >= 0) & (x <= 1) & (z >= 0) & (z <= 1)
    sx, sz = np.sin(PI * x), np.sin(PI * z)
    return (np.where(inside, -sx * sx * np.sin(2 * PI * z), 0.0),
            np.where(inside, np.sin(2 * PI * x) * sz * sz, 0.0))


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


def wmax1_2d(T, period, revolutions, amp, centre=(0.5, 0.5), n=900, nt=1200):
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
        best = max(best, float(np.max(np.abs(Ux) + np.abs(Uz))))
    return best


def wmax1_3d(T, period, revolutions, field, centre=(0.5, 0.5), swirl=(0.5, 0.5),
             n=170, nt=420):
    g = (np.arange(n) + 0.5) / n
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    thd = 2 * PI * revolutions / T
    best = 0.0
    for t in np.linspace(0.0, T, nt):
        th = 2 * PI * revolutions * t / T
        ct, st = math.cos(th), math.sin(th)
        a = math.cos(2 * PI * t / period) if period > 0 else 1.0
        rx, ry = X - centre[0], Y - centre[1]
        xix, xiy = centre[0] + (rx * ct + ry * st), centre[1] + (-rx * st + ry * ct)
        ux, uy, uz = u0_3d(xix, xiy, Z, field, swirl)
        ux, uy, uz = a * ux, a * uy, a * uz
        Ux = thd * (-ry) + (ux * ct - uy * st)
        Uy = thd * (rx) + (ux * st + uy * ct)
        best = max(best, float(np.max(np.abs(Ux) + np.abs(Uy) + np.abs(uz))))
    return best


def steps_for(N, co, u_max, end_time, multiple=4):
    """The step count build_case writes for a convergence case."""
    dt_target = co / (float(N) * u_max)
    return multiple * int(math.ceil(end_time / (multiple * dt_target)))


def main():
    import argparse
    import json

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--field", required=True, choices=["vortex", "spiralling", "leveque"])
    ap.add_argument("--end-time", type=float, required=True)
    ap.add_argument("--period", type=float, required=True)
    ap.add_argument("--revolutions", type=float, required=True)
    ap.add_argument("--translation", type=float, default=0.0)
    ap.add_argument("--centre", required=True, help='rotation centre, e.g. "(0.5 0 0.5)"')
    ap.add_argument("--swirl", default="(0.5 0.5 0)", help="3D shear swirl axis")
    ap.add_argument("--max-alpha-co", type=float, required=True)
    ap.add_argument("--resolutions", required=True, help="comma-separated N")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    c = [float(v) for v in a.centre.strip("()").split()]
    sw = [float(v) for v in a.swirl.strip("()").split()]
    if a.field == "vortex":
        centre = (c[0], c[2])                  # the 2D case lives in the x-z plane
        framed = wmax1_2d(a.end_time, a.period, a.revolutions, a.translation, centre)
        plain = wmax1_2d(a.end_time, a.period, 0.0, 0.0, centre)
    else:
        if a.translation != 0.0:
            raise SystemExit("the 3D bound has no translating frame")
        centre = (c[0], c[1])
        framed = wmax1_3d(a.end_time, a.period, a.revolutions, a.field, centre, sw[:2])
        plain = wmax1_3d(a.end_time, a.period, 0.0, a.field, centre, sw[:2])

    Ns = [int(n) for n in a.resolutions.split(",")]
    out = {
        "u_max": framed,
        "u_max_unframed": plain,
        "norm": "max over domain and period of ||w||_1",
        "field": a.field,
        "end_time": a.end_time,
        "period": a.period,
        "revolutions": a.revolutions,
        "translation": a.translation,
        "rotation_centre": list(centre),
        "max_alpha_co": a.max_alpha_co,
        "convergence_steps": {str(N): steps_for(N, a.max_alpha_co, framed, a.end_time)
                              for N in Ns},
    }
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
