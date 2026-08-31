#!/usr/bin/env python3
"""3D interface snapshots at t = 0, T/4, T/2, 3T/4, T.

The volume fraction is read with foamlib, reshaped onto the uniform hex block
and contoured at alpha = 0.5 with marching cubes. The framed case is drawn
twice: in the laboratory frame, and mapped back through E_t^{-1} into the
co-rotating frame, where it must reproduce the unframed case.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import re

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
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure

FACE = "#3d6fb4"
EDGE = "#12263f"

# domain [0,1] x [0,1] x [0,LZ]; set from the command line
LX = LY = 1.0
LZ = 2.0
# the drawn box may be shorter than the domain when the material never
# reaches the top; set from the command line, 0 means 'no clip'
ZCLIP = 0.0
NZF = 2


def read_alpha(case_dir: str, t: float, N: int) -> np.ndarray:
    """Volume fraction as a (nz, ny, nx) array; blockMesh orders x fastest."""
    tn = "%g" % t
    for name in ("alpha.water", "alpha.water.gz"):
        p = os.path.join(case_dir, tn, name)
        if os.path.exists(p):
            op = gzip.open if p.endswith(".gz") else open
            with op(p, "rt") as f:
                s = f.read()
            break
    else:
        raise SystemExit(f"no alpha.water at t={tn} in {case_dir}")

    i = s.index("internalField")
    if "nonuniform" not in s[i:i + 60]:
        val = float(s[i:s.index(";", i)].split()[-1])
        return np.full((NZF * N, N, N), val)
    j = s.index("(", i)
    k = s.index(")", j)
    a = np.fromstring(s[j + 1:k], sep=" ")
    nx = ny = N
    nz = NZF * N
    if a.size != nx * ny * nz:
        raise SystemExit(f"expected {nx*ny*nz} values, got {a.size}")
    return a.reshape(nz, ny, nx)


def isosurface(vol: np.ndarray, N: int):
    """alpha = 0.5 surface as (verts_xyz, faces); empty if the phase vanished."""
    nz, ny, nx = vol.shape
    if vol.min() > 0.5 or vol.max() < 0.5:
        return None, None
    # pad so surfaces touching the block boundary still close
    padded = np.pad(vol, 1, mode="constant", constant_values=0.0)
    dz, dy, dx = LZ / nz, LY / ny, LX / nx
    verts, faces, _, _ = measure.marching_cubes(
        padded, level=0.5, spacing=(dz, dy, dx)
    )
    # undo the pad offset and the half-cell shift of the cell-centred data
    verts = verts - np.array([dz, dy, dx])
    verts = verts + np.array([0.5 * dz, 0.5 * dy, 0.5 * dx])
    # (z, y, x) -> (x, y, z)
    return verts[:, ::-1], faces


def to_moving_frame(verts, t, T, centre, revolutions):
    """E_t^{-1}(x) = c + Q(t)^T (x - c), rotation about the axis along e_z."""
    th = 2.0 * math.pi * revolutions * t / T
    ct, st = math.cos(th), math.sin(th)
    dx = verts[:, 0] - centre[0]
    dy = verts[:, 1] - centre[1]
    out = verts.copy()
    out[:, 0] = centre[0] + ct * dx + st * dy
    out[:, 1] = centre[1] - st * dx + ct * dy
    return out


def _nz_label(nzf):
    """'N' rather than '1N' when the cells are already cubic."""
    return "N" if nzf == 1 else f"{nzf}N"

def draw(ax, verts, faces, title, centre=None):
    if verts is not None and len(faces):
        mesh = Poly3DCollection(verts[faces], alpha=0.9)
        mesh.set_facecolor(FACE)
        mesh.set_edgecolor(EDGE)
        mesh.set_linewidth(0.05)
        ax.add_collection3d(mesh)
    if centre is not None:
        ax.plot([centre[0]], [centre[1]], [0.0], marker="+", color="k", ms=12)
    ax.set_xlim(0, LX)
    ax.set_ylim(0, LY)
    zt = ZCLIP if ZCLIP > 0 else LZ
    ax.set_zlim(0, zt)
    ax.set_box_aspect((1, 1, zt))
    # x's last tick and y's first tick meet at the corner nearest the viewer,
    # and both read the same two digits in every panel; the axis names carry
    # more and cannot collide.
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["", ""])
    ax.set_yticklabels(["", ""])
    ax.set_xlabel("$x$", labelpad=-8)
    ax.set_ylabel("$y$", labelpad=-8)
    if zt > 1.5:
        ax.set_zticks([0, 1, 2])
    elif zt >= 1.0:
        ax.set_zticks([0, 0.5, 1])
    else:
        ax.set_zticks([0, round(0.5 * zt, 2), zt])
    ax.tick_params(axis='z', labelsize=FS - 6, pad=2)
    ax.tick_params(axis='x', labelsize=FS - 6)
    ax.tick_params(axis='y', labelsize=FS - 6)
    ax.view_init(elev=18, azim=-58)
    if title:
        ax.set_title(title, fontsize=FS, pad=0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--framed-case", required=True)
    p.add_argument("--plain-case", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--scheme", required=True)
    p.add_argument("--end-time", type=float, required=True)
    p.add_argument("--centre", required=True)
    p.add_argument("--revolutions", type=float, required=True)
    p.add_argument("--lz", type=float, default=2.0)
    p.add_argument("--nz-factor", type=int, default=2)
    p.add_argument("--z-clip", type=float, default=0.0,
                   help="draw the box only up to this height (0 = full)")
    p.add_argument("--field", default="spiralling deformation")
    p.add_argument("--out-moving", required=True)
    p.add_argument("--out-static", required=True)
    a = p.parse_args()

    global LZ, NZF, ZCLIP
    LZ = a.lz
    ZCLIP = a.z_clip
    NZF = a.nz_factor
    T, N = a.end_time, a.n
    c = [float(v) for v in a.centre.strip("()").split()]
    centre = (c[0], c[1])
    times = [0.0, T / 4, T / 2, 3 * T / 4, T]
    labels = ["$t=0$", "$t=T/4$", "$t=T/2$", "$t=3T/4$", "$t=T$"]

    for out in (a.out_moving, a.out_static):
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)

    # ---- framed case, laboratory and co-rotating views ----
    fig = plt.figure(figsize=(30.0, 15.0))
    for col, (t, lab) in enumerate(zip(times, labels)):
        v, f = isosurface(read_alpha(a.framed_case, t, N), N)
        ax = fig.add_subplot(2, len(times), col + 1, projection="3d")
        draw(ax, v, f, lab, centre=centre)
        if col == 0:
            ax.text2D(-0.16, 0.5, "laboratory frame\n(as simulated)",
                      transform=ax.transAxes, rotation=90, va="center", fontsize=FS - 2)
        ax2 = fig.add_subplot(2, len(times), len(times) + col + 1, projection="3d")
        draw(ax2,
             None if v is None else to_moving_frame(v, t, T, centre, a.revolutions),
             f, "")
        if col == 0:
            ax2.text2D(-0.16, 0.5, "co-rotating frame\n$E_t^{-1}$ applied",
                       transform=ax2.transAxes, rotation=90, va="center", fontsize=FS - 2)
    fig.suptitle(
        f"Moving-frame case: {a.revolutions:g} revolution about the axis through "
        f"$({centre[0]:g},{centre[1]:g})$ superposed on the {a.field}\n"
        f"{a.scheme}, $N={N}$ ($n_z={_nz_label(NZF)}$), $T={T:g}$",
        fontsize=FS + 2,
    )
    fig.tight_layout(rect=[0.02, 0.0, 1, 0.93])
    fig.savefig(a.out_moving, dpi=170)
    plt.close(fig)
    print(a.out_moving)

    # ---- unframed case ----
    fig = plt.figure(figsize=(30.0, 8.0))
    for col, (t, lab) in enumerate(zip(times, labels)):
        v, f = isosurface(read_alpha(a.plain_case, t, N), N)
        ax = fig.add_subplot(1, len(times), col + 1, projection="3d")
        draw(ax, v, f, lab)
        if col == 0:
            ax.text2D(-0.16, 0.5, "laboratory frame\n(= co-rotating)",
                      transform=ax.transAxes, rotation=90, va="center", fontsize=FS - 2)
    fig.suptitle(
        f"Non-moving-frame case: the original {a.field} "
        "(no frame superposed)\n"
        f"{a.scheme}, $N={N}$; the motion for $t>T/2$ retraces the motion for $t<T/2$",
        fontsize=FS + 2,
    )
    fig.tight_layout(rect=[0.02, 0.0, 1, 0.88])
    fig.savefig(a.out_static, dpi=170)
    plt.close(fig)
    print(a.out_static)


if __name__ == "__main__":
    main()
