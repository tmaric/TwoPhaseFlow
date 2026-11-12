#!/usr/bin/env python3
"""Interface snapshots at t = 0, T/2, T, read with foamlib.

The framed case is drawn twice: in the laboratory frame, i.e. as simulated, and
mapped back through E_t^{-1} into the co-rotating frame, where it must reproduce
the unframed case. That agreement is the frame-consistency check.

Three snapshots, not five, and no title: the figures are printed at \\textwidth,
and the review asked for larger panels and fonts. The caption carries what the
title used to say.
"""

from __future__ import annotations

import argparse
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Reviewer comment: figures are reduced heavily in print, so the default
# matplotlib sizes are unreadable. Roughly double them, and grow the canvas to
# match so labels do not collide.
FS = 20
plt.rcParams.update({
    "font.size": FS,
    "axes.titlesize": FS + 2,
    "axes.labelsize": FS,
    "xtick.labelsize": FS - 4,
    "ytick.labelsize": FS - 4,
    "legend.fontsize": FS - 4,
})

from foamlib import FoamCase

FILL = "#3d6fb4"
EDGE = "#12263f"
REF = "#c0392b"

DISC_CENTRE = (0.5, 0.75)
DISC_RADIUS = 0.15


def alpha_at(case_dir: str, t: float, N: int) -> np.ndarray:
    """Volume fraction as an (N, N) array indexed [z, x]."""
    case = FoamCase(case_dir)
    chosen = None
    for td in case:
        if abs(float(td.time) - t) < 1e-6:
            chosen = td
            break
    if chosen is None:
        raise SystemExit(
            f"time {t} not written in {case_dir}; have "
            + ", ".join(str(td.time) for td in case)
        )
    field = chosen["alpha.water"].internal_field
    a = np.asarray(field, dtype=float).ravel()
    if a.size == 1:                     # uniform
        a = np.full(N * N, a[0])
    if a.size != N * N:
        raise SystemExit(f"expected {N*N} values, got {a.size} in {case_dir} t={t}")
    # blockMesh orders a single hex block with x fastest, then y, then z
    return a.reshape(N, N)


def contours(a: np.ndarray, N: int, level: float = 0.5):
    x = (np.arange(N) + 0.5) / N
    X, Z = np.meshgrid(x, x)
    fig = plt.figure()
    cs = plt.contour(X, Z, a, levels=[level])
    segs = [np.asarray(s) for s in cs.allsegs[0]]
    plt.close(fig)
    return segs


def to_moving_frame(seg, t, T, centre, revolutions):
    """E_t^{-1}(x) = c + Q(t)^T (x - c)."""
    th = 2.0 * math.pi * revolutions * t / T
    ct, st = math.cos(th), math.sin(th)
    dx = seg[:, 0] - centre[0]
    dz = seg[:, 1] - centre[1]
    return np.column_stack(
        (centre[0] + ct * dx + st * dz, centre[1] - st * dx + ct * dz)
    )


def draw(ax, segs, title, centre=None):
    for s in segs:
        ax.fill(s[:, 0], s[:, 1], color=FILL, alpha=0.85, zorder=2)
        ax.plot(s[:, 0], s[:, 1], color=EDGE, lw=1.6, zorder=3)
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(
        DISC_CENTRE[0] + DISC_RADIUS * np.cos(th),
        DISC_CENTRE[1] + DISC_RADIUS * np.sin(th),
        color=REF, lw=2.0, ls="--", zorder=4,
    )
    if centre is not None:
        ax.plot(*centre, marker="+", color="k", ms=14, mew=2.4, zorder=5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
    ax.set_yticks([0, 0.5, 1], ["0", "0.5", "1"])
    ax.tick_params(labelsize=FS - 4, length=2)
    ax.set_title(title, fontsize=FS)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)


def inner_ticks_off(axes):
    """Tick labels on the outer panels only; the panels share their axes."""
    nr, nc = axes.shape
    for i in range(nr):
        for j in range(nc):
            if j > 0:
                axes[i, j].set_yticklabels([])
            if i < nr - 1:
                axes[i, j].set_xticklabels([])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--framed-case", required=True)
    p.add_argument("--plain-case", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--scheme", required=True)
    p.add_argument("--end-time", type=float, required=True)
    p.add_argument("--centre", required=True, help='e.g. "(0.5 0 0.5)"')
    p.add_argument("--revolutions", type=float, required=True)
    p.add_argument("--out-moving", required=True)
    p.add_argument("--out-static", required=True)
    args = p.parse_args()

    T = args.end_time
    N = args.n
    c = [float(v) for v in args.centre.strip("()").split()]
    centre = (c[0], c[2])
    times = [0.0, T / 2, T]
    labels = ["$t=0$", "$t=T/2$", "$t=T$"]

    for out in (args.out_moving, args.out_static):
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)

    # ---- framed case: laboratory frame and co-rotating frame ----
    # 12 in wide at \textwidth (about 6 in) halves every size: FS = 20 prints
    # at 10 pt, the caption size
    fig, axes = plt.subplots(2, len(times), figsize=(12.0, 8.4))
    for col, (t, lab) in enumerate(zip(times, labels)):
        segs = contours(alpha_at(args.framed_case, t, N), N)
        draw(axes[0, col], segs, lab, centre=centre)
        draw(
            axes[1, col],
            [to_moving_frame(s, t, T, centre, args.revolutions) for s in segs],
            "",
        )
    axes[0, 0].set_ylabel("laboratory frame\n(as simulated)", fontsize=FS)
    axes[1, 0].set_ylabel("rotating frame\n(mapped back)", fontsize=FS)
    inner_ticks_off(axes)
    fig.tight_layout(h_pad=0.6, w_pad=0.4)
    fig.savefig(args.out_moving, dpi=200)
    plt.close(fig)
    print(args.out_moving)

    # ---- unframed case ----
    fig, axes = plt.subplots(1, len(times), figsize=(12.0, 4.5))
    for col, (t, lab) in enumerate(zip(times, labels)):
        draw(axes[col], contours(alpha_at(args.plain_case, t, N), N), lab)
    axes[0].set_ylabel("laboratory frame", fontsize=FS)
    inner_ticks_off(axes[None, :])
    fig.tight_layout(w_pad=0.4)
    fig.savefig(args.out_static, dpi=200)
    plt.close(fig)
    print(args.out_static)


if __name__ == "__main__":
    main()
