#!/usr/bin/env python3
"""Interface evolution every 0.1T, with and without the frame, from simulation data.

This is the marker figure of the Bothe & Maric recipe redrawn from the solver's
output. Rows are the original reversed vortex and the same base field with the
frame superposed; columns are the two half-periods. A snapshot at t and its
mirror at T-t share a colour, so retracing is read along a row: in the first row
the right panel repeats the left one curve for curve, in the second it does not.

Simulated interfaces are the alpha = 0.5 contours. The exact interface under
each is computed here, independently of the solver: the base flow is
alpha(t) u0, so a base-frame trajectory is the autonomous flow of u0 run to the
pseudo-time s(t) = int_0^t alpha = (P/2pi) sin(2pi t/P), and the laboratory
position is c + Q(t)(y - c). Where a coloured curve leaves its dark line, that
is solver error.

The same two runs are also written as an animation, side by side.
"""

from __future__ import annotations

import argparse
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

from visualize import DISC_CENTRE, DISC_RADIUS, alpha_at, contours

# Same reviewer-driven sizes as the other snapshot figures.
FS = 20
plt.rcParams.update({
    "font.size": FS,
    "axes.titlesize": FS + 2,
    "axes.labelsize": FS,
    "xtick.labelsize": FS - 6,
    "ytick.labelsize": FS - 6,
    "legend.fontsize": FS - 4,
})

INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#d9d8d3"
# Pair index j = min(k, 10-k) for t = 0.1kT. j = 0 is the disc itself (t = 0
# and t = T), drawn in ink; j = 1..5 take the first five categorical slots in
# fixed order. Validated (light): CVD and normal-vision separation pass; three
# slots sit below 3:1 on the surface, so every curve is also labelled directly.
PAIR_COLOURS = [INK, "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
ALPHA_CMAP = LinearSegmentedColormap.from_list("alpha", ["#fcfcfb", "#256abf"])


# ------------------------------------------------------------------ exact
def u0(p: np.ndarray) -> np.ndarray:
    """The reversed-vortex base field in the (x, z) plane, zero outside [0,1]^2."""
    x, z = p[:, 0], p[:, 1]
    inside = (x >= 0) & (x <= 1) & (z >= 0) & (z <= 1)
    ux = -np.sin(np.pi * x) ** 2 * np.sin(2 * np.pi * z)
    uz = np.sin(2 * np.pi * x) * np.sin(np.pi * z) ** 2
    return np.where(inside[:, None], np.stack([ux, uz], 1), 0.0)


def exact_interfaces(times, T, period, centre, revolutions, n=6000, ds=2e-3):
    """Exact laboratory-frame boundary of the disc at each requested time."""
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    y = np.stack([DISC_CENTRE[0] + DISC_RADIUS * np.cos(a),
                  DISC_CENTRE[1] + DISC_RADIUS * np.sin(a)], 1)
    s_of = {t: period / (2 * np.pi) * math.sin(2 * np.pi * t / period) for t in times}
    base = {}
    s_now = 0.0
    for s_target in sorted(set(s_of.values())):
        steps = max(1, int(math.ceil((s_target - s_now) / ds)))
        h = (s_target - s_now) / steps
        for _ in range(steps):
            k1 = u0(y)
            k2 = u0(y + 0.5 * h * k1)
            k3 = u0(y + 0.5 * h * k2)
            k4 = u0(y + h * k3)
            y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        s_now = s_target
        base[s_target] = y.copy()
    out = {}
    for t in times:
        yb = base[s_of[t]]
        th = 2 * np.pi * revolutions * t / T
        ct, st = math.cos(th), math.sin(th)
        r = yb - centre
        lab = np.column_stack((centre[0] + ct * r[:, 0] - st * r[:, 1],
                               centre[1] + st * r[:, 0] + ct * r[:, 1]))
        out[t] = np.vstack([lab, lab[:1]])
    return out


# ------------------------------------------------------------------ drawing
def frame_axes(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])
    ax.tick_params(length=2, colors=MUTED)
    ax.grid(color=GRID, lw=0.6, zorder=0)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)
        sp.set_color(MUTED)


def label_points(curves):
    """For each curve, the point farthest from every other curve.

    Labels go where a curve is most alone, so they rarely sit on a neighbour,
    and away from labels already placed, so two labels do not overlap.
    """
    pts = [np.vstack(c)[:: max(1, len(np.vstack(c)) // 400)] for c in curves]
    out = []
    for i, p in enumerate(pts):
        others = [q for j, q in enumerate(pts) if j != i]
        if not others:
            out.append(p[len(p) // 2])
            continue
        o = np.vstack(others)
        d = np.min(np.hypot(p[:, None, 0] - o[None, :, 0],
                            p[:, None, 1] - o[None, :, 1]), axis=1)
        if out:
            # a label is about 0.08 wide; closer than that to a placed one overlaps
            q = np.array(out)
            dl = np.min(np.hypot(p[:, None, 0] - q[None, :, 0],
                                 p[:, None, 1] - q[None, :, 1]), axis=1)
            d = np.where(dl > 0.1, d, d - 1.0)
        # keep labels off the domain edge
        d = np.where((p[:, 0] > 0.06) & (p[:, 0] < 0.94) &
                     (p[:, 1] > 0.06) & (p[:, 1] < 0.94), d, -2.0)
        out.append(p[int(np.argmax(d))])
    return out


def panel(ax, ks, sim, exact, centre=None):
    frame_axes(ax)
    curves = []
    for k in ks:
        j = min(k, 10 - k)
        col = PAIR_COLOURS[j]
        for seg in exact[k]:
            ax.plot(seg[:, 0], seg[:, 1], color=MUTED, lw=0.8, zorder=2)
        for seg in sim[k]:
            ax.plot(seg[:, 0], seg[:, 1], color=col, lw=2.2, zorder=3)
        curves.append(sim[k])
    for k, (x, y) in zip(ks, label_points(curves)):
        txt = "$0$" if k == 0 else ("$T$" if k == 10 else f"${k / 10:g}T$")
        ax.annotate(txt, (x, y), xytext=(4, 4), textcoords="offset points",
                    fontsize=FS - 6, color=INK, zorder=5,
                    bbox=dict(boxstyle="round,pad=0.15", fc="#fcfcfb", ec="none", alpha=0.85))
    if centre is not None:
        ax.plot(*centre, marker="+", color=INK, ms=14, mew=2.0, zorder=6)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--framed-case", required=True)
    p.add_argument("--plain-case", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--scheme", required=True)
    p.add_argument("--end-time", type=float, required=True)
    p.add_argument("--period", type=float, required=True)
    p.add_argument("--centre", required=True, help='e.g. "(0.5 0 0.5)"')
    p.add_argument("--revolutions", type=float, required=True)
    p.add_argument("--out", required=True, help="PNG; a PDF is written beside it")
    p.add_argument("--out-anim", required=True, help="animated GIF")
    p.add_argument("--anim-interval", type=float, default=None,
                   help="time between animation frames (default: every written time)")
    a = p.parse_args()

    T, N = a.end_time, a.n
    c3 = [float(v) for v in a.centre.strip("()").split()]
    centre = np.array([c3[0], c3[2]])
    snap = [k * T / 10 for k in range(11)]

    for out in (a.out, a.out_anim):
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)

    # ---- the static figure: snapshots every 0.1T ----
    ex_plain = exact_interfaces(snap, T, a.period, centre, 0.0)
    ex_framed = exact_interfaces(snap, T, a.period, centre, a.revolutions)
    sim_plain = {k: contours(alpha_at(a.plain_case, t, N), N) for k, t in enumerate(snap)}
    sim_framed = {k: contours(alpha_at(a.framed_case, t, N), N) for k, t in enumerate(snap)}
    ex_plain = {k: [ex_plain[t]] for k, t in enumerate(snap)}
    ex_framed = {k: [ex_framed[t]] for k, t in enumerate(snap)}

    first, second = list(range(0, 6)), list(range(5, 11))
    fig, axes = plt.subplots(2, 2, figsize=(15.5, 16.0))
    panel(axes[0, 0], first, sim_plain, ex_plain)
    panel(axes[0, 1], second, sim_plain, ex_plain)
    panel(axes[1, 0], first, sim_framed, ex_framed, centre)
    panel(axes[1, 1], second, sim_framed, ex_framed, centre)
    axes[0, 0].set_title(r"$0 \leq t \leq T/2$")
    axes[0, 1].set_title(r"$T/2 \leq t \leq T$")
    axes[0, 0].set_ylabel("original\n(retracing)")
    axes[1, 0].set_ylabel(f"frame superposed\n({a.revolutions:g} revolution about "
                          f"$({centre[0]:g},{centre[1]:g})$)")
    for ax in axes[1]:
        ax.set_xlabel("$x$")

    handles = [Line2D([], [], color=PAIR_COLOURS[j], lw=3,
                      label=("$t=0,\\ T$" if j == 0 else f"$t={j / 10:g}T,\\ {1 - j / 10:g}T$"
                             if j < 5 else "$t=0.5T$"))
               for j in range(6)]
    handles.append(Line2D([], [], color=MUTED, lw=0.8, label="exact"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(f"Reversed vortex, {a.scheme}, $N={N}$: simulated $\\alpha=0.5$ "
                 "every $0.1T$ (colour), exact interface (thin)", fontsize=FS)
    fig.tight_layout(rect=[0, 0.07, 1, 0.97])
    fig.savefig(a.out, dpi=160)
    # contours are line art: a vector copy for the manuscript
    pdf = os.path.splitext(a.out)[0] + ".pdf"
    fig.savefig(pdf)
    plt.close(fig)
    print(a.out)
    print(pdf)

    # ---- the animation: every written time, both runs side by side ----
    from foamlib import FoamCase

    def written(case):
        return {round(float(td.time), 9) for td in FoamCase(case)}

    times = sorted(written(a.plain_case) & written(a.framed_case))
    if a.anim_interval:
        step = a.anim_interval
        times = [t for t in times if abs(t / step - round(t / step)) < 1e-6]
    ex_p = exact_interfaces(times, T, a.period, centre, 0.0)
    ex_f = exact_interfaces(times, T, a.period, centre, a.revolutions)
    fields = {t: (alpha_at(a.plain_case, t, N), alpha_at(a.framed_case, t, N)) for t in times}

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 7.2))
    ext = (0, 1, 0, 1)

    def draw(i):
        t = times[i]
        ap, af = fields[t]
        for ax, alpha, ex, title, cen in (
                (axes[0], ap, ex_p[t], "original (retracing)", None),
                (axes[1], af, ex_f[t], "frame superposed", centre)):
            ax.clear()
            frame_axes(ax)
            ax.imshow(alpha, origin="lower", extent=ext, cmap=ALPHA_CMAP,
                      vmin=0, vmax=1, interpolation="nearest", zorder=1)
            ax.plot(ex[:, 0], ex[:, 1], color=INK, lw=0.9, zorder=3)
            ax.set_title(title, fontsize=FS - 2)
            if cen is not None:
                th = 2 * np.pi * a.revolutions * t / T
                r = 0.47
                # a spoke turning with the frame, so the rotation is visible
                ax.plot([cen[0], cen[0] + r * math.cos(th)],
                        [cen[1], cen[1] + r * math.sin(th)],
                        color=MUTED, lw=1.2, ls="--", zorder=2)
                ax.plot(*cen, marker="+", color=INK, ms=12, mew=2.0, zorder=4)
        fig.suptitle(f"$t = {t / T:.3f}\\,T$    ({a.scheme}, $N={N}$; "
                     "shading: simulated $\\alpha$, line: exact)", fontsize=FS - 2)
        return []

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    anim = FuncAnimation(fig, draw, frames=len(times), blit=False)
    anim.save(a.out_anim, writer=PillowWriter(fps=6), dpi=80)
    plt.close(fig)
    print(a.out_anim)


if __name__ == "__main__":
    main()
