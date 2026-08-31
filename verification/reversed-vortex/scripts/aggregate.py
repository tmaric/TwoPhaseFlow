"""Collect the final-time errors of every case into a CSV and a LaTeX table.

Convergence order is computed for the L1 geometrical shape error only; the
volume and boundedness errors are tabulated without an order, because they do
not converge in the usual sense -- they sit at round-off.
"""

from __future__ import annotations

import csv
import math
import os

FRAME_LABEL = {"none": "no frame (retracing)", "frameA": "Frame A (non-retracing)"}


def _read_final(path: str, end_time: float):
    """Return (Eshape, Emass, Ebound) at t == end_time, or None."""
    if not os.path.exists(path):
        return None
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        r = line.split()
        if abs(float(r[0]) - end_time) < 1e-9:
            return float(r[1]), float(r[2]), float(r[3])
    return None


def _orders(ns, es):
    out = [None]
    for i in range(1, len(ns)):
        if es[i - 1] and es[i]:
            out.append(math.log(es[i - 1] / es[i]) / math.log(ns[i] / ns[i - 1]))
        else:
            out.append(None)
    return out


def _fmt(x):
    """1.234e-03 -> $1.234\\times10^{-3}$"""
    if x is None:
        return "---"
    if x == 0.0:
        return "$0$"
    m, e = f"{x:.3e}".split("e")
    return rf"${m}\times10^{{{int(e)}}}$"


def write(runs, schemes, frames, resolutions, end_time, csv_path, tex_path,
          test_name="reversed-vortex", label_prefix="conv", workflow="reversed-vortex"):
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)

    rows = []
    for scheme in schemes:
        for frame in frames:
            for N in resolutions:
                case = f"{scheme}_{frame}_N{N}"
                v = _read_final(
                    os.path.join(runs, case, "volumeFractionError.dat"), end_time
                )
                if v is None:
                    continue
                rows.append(
                    {
                        "scheme": scheme,
                        "frame": frame,
                        "N": N,
                        "E_shape": v[0],
                        "E_mass": v[1],
                        "E_bound": v[2],
                    }
                )

    # observed orders, per (scheme, frame) series
    for scheme in schemes:
        for frame in frames:
            sel = sorted(
                [r for r in rows if r["scheme"] == scheme and r["frame"] == frame],
                key=lambda r: r["N"],
            )
            ns = [r["N"] for r in sel]
            es = [r["E_shape"] for r in sel]
            for r, o in zip(sel, _orders(ns, es)):
                r["order"] = o

    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["scheme", "frame", "N", "E_shape", "order", "E_mass", "E_bound"]
        )
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in w.fieldnames})

    # ---- LaTeX ----------------------------------------------------------
    L = []
    for scheme in schemes:
        by = {
            frame: {r["N"]: r for r in rows if r["scheme"] == scheme and r["frame"] == frame}
            for frame in frames
        }
        present = sorted({N for frame in frames for N in by[frame]})
        if not present:
            continue
        L.append(r"\begin{table}[htbp]")
        L.append(r"\centering")
        L.append(r"\small")
        L.append(r"\setlength{\tabcolsep}{4pt}")
        L.append(
            rf"\caption{{{scheme} reconstruction on the {test_name} test. "
            r"The convergence order is reported for the L1 geometrical shape "
            r"error only. Data: "
            rf"\cite{{figshare2026}}, \path{{{workflow}/results/convergence.csv}}.}}"
        )
        L.append(rf"\label{{tab:{label_prefix}-{scheme}}}")
        L.append(r"\begin{tabular}{r cc c cc c}")
        L.append(r"\toprule")
        L.append(
            r"& \multicolumn{3}{c}{"
            + FRAME_LABEL[frames[0]]
            + r"} & \multicolumn{3}{c}{"
            + FRAME_LABEL[frames[1]]
            + r"} \\"
        )
        L.append(r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}")
        L.append(
            r"$N$ & $E_{\text{shape}}$ & order & $E_{\text{bound}}$"
            r" & $E_{\text{shape}}$ & order & $E_{\text{bound}}$ \\"
        )
        L.append(r"\midrule")
        for N in present:
            cells = [str(N)]
            for frame in frames:
                r = by[frame].get(N)
                if r is None:
                    cells += ["---", "---", "---"]
                else:
                    o = r.get("order")
                    cells += [
                        _fmt(r["E_shape"]),
                        "---" if o is None else f"{o:.2f}",
                        _fmt(r["E_bound"]),
                    ]
            L.append(" & ".join(cells) + r" \\")
        L.append(r"\bottomrule")
        L.append(r"\end{tabular}")
        L.append(r"\end{table}")
        L.append("")

    with open(tex_path, "w") as f:
        f.write("\n".join(L))
