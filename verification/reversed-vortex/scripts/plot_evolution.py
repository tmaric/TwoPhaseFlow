#!/usr/bin/env python3
"""Retracing, seen in the interface: one mirror pair, with and without the frame.

The static figure says one thing. Without the frame, the interface at T - t lies
on the interface at t, because the second half-period retraces the first; with
the frame, the two are apart. Each panel overlays the simulated alpha = 0.5
contours at 0.3T (solid) and its mirror time 0.7T (dashed). An earlier version
drew every 0.1T against the exact interface; the review found it overloaded.

The exact interface used by the animation is computed here, independently of
the solver.

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
# Categorical slots validated for CVD and normal-vision separation (light). The
# mirror pair uses slots 1 and 2, and also differs in line width and dash, so
# it does not rely on colour alone.
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
    p.add_argument("--mirror", type=float, default=0.3,
                   help="draw the pair t = mirror*T and (1 - mirror)*T")
    p.add_argument("--static-only", action="store_true",
                   help="skip the animation")
    a = p.parse_args()

    T, N = a.end_time, a.n
    c3 = [float(v) for v in a.centre.strip("()").split()]
    centre = np.array([c3[0], c3[2]])

    for out in (a.out, a.out_anim):
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)

    # ---- the static figure: one mirror pair per run, overlaid ----
    m = a.mirror
    t1, t2 = m * T, (1.0 - m) * T
    # 12 in wide at \textwidth halves every size: FS = 20 prints at 10 pt
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 6.6))
    for ax, case, title, cen in (
            (axes[0], a.plain_case, "without frame", None),
            (axes[1], a.framed_case, "with frame", centre)):
        frame_axes(ax)
        ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
        ax.set_yticks([0, 0.5, 1], ["0", "0.5", "1"])
        ax.tick_params(labelsize=FS - 4)
        for seg in contours(alpha_at(case, t1, N), N):
            ax.plot(seg[:, 0], seg[:, 1], color=PAIR_COLOURS[1], lw=5.0,
                    solid_capstyle="round", zorder=3)
        for seg in contours(alpha_at(case, t2, N), N):
            ax.plot(seg[:, 0], seg[:, 1], color=PAIR_COLOURS[2], lw=2.2,
                    ls=(0, (4, 2.5)), zorder=4)
        if cen is not None:
            ax.plot(*cen, marker="+", color=INK, ms=16, mew=2.4, zorder=6)
        ax.set_title(title, fontsize=FS + 2)
    axes[1].set_yticklabels([])
    handles = [Line2D([], [], color=PAIR_COLOURS[1], lw=5.0, label=f"$t={m:g}\\,T$"),
               Line2D([], [], color=PAIR_COLOURS[2], lw=2.2, ls=(0, (4, 2.5)),
                      label=f"$t={1 - m:g}\\,T$")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=FS, bbox_to_anchor=(0.5, 0.0), handlelength=3.0)
    fig.tight_layout(rect=[0, 0.08, 1, 1], w_pad=1.0)
    fig.savefig(a.out, dpi=160)
    # contours are line art: a vector copy for the manuscript
    pdf = os.path.splitext(a.out)[0] + ".pdf"
    fig.savefig(pdf)
    plt.close(fig)
    print(a.out)
    print(pdf)
    if a.static_only:
        return

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
