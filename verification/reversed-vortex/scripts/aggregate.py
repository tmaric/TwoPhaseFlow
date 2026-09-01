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


def _fit_order(ns, es):
    """Least-squares order over the whole sequence, and its residual.

    Returns (order, max |residual| in log E, R^2). The residual is reported
    because a sequence can have a respectable fitted order while departing
    badly from a power law, and that departure is itself a result.
    """
    pts = [(n, e) for n, e in zip(ns, es) if n and e and e > 0]
    if len(pts) < 3:
        return None, None, None
    xs = [math.log(float(n)) for n, _ in pts]
    ys = [math.log(e) for _, e in pts]
    n = len(pts)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0.0:
        return None, None, None
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    icpt = my - slope * mx
    res = [y - (slope * x + icpt) for x, y in zip(xs, ys)]
    ss_res = sum(r * r for r in res)
    ss_tot = sum((y - my) ** 2 for y in ys)
    r2 = 1.0 - ss_res / ss_tot if ss_tot else None
    return -slope, max(abs(r) for r in res), r2


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
                    # a level silently missing from the table is worse than a
                    # failed build: it looks like a study that was never run
                    raise SystemExit(
                        f"aggregate: {case} has no row at t={end_time} in "
                        f"volumeFractionError.dat -- the run did not finish, or "
                        f"its output was not flushed. Refusing to write a table "
                        f"with a missing level."
                    )
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

    # the fit is a property of a whole series, so it repeats down each one
    fits = {}
    for scheme in schemes:
        for frame in frames:
            sel = sorted(
                [r for r in rows if r["scheme"] == scheme and r["frame"] == frame],
                key=lambda r: r["N"],
            )
            fits[(scheme, frame)] = _fit_order(
                [r["N"] for r in sel], [r["E_shape"] for r in sel]
            )
    for r in rows:
        o, res, r2 = fits[(r["scheme"], r["frame"])]
        r["fit_order"] = o
        r["fit_residual"] = res
        r["fit_r2"] = r2

    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["scheme", "frame", "N", "E_shape", "order", "E_mass",
                           "E_bound", "fit_order", "fit_residual", "fit_r2"]
        )
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in w.fieldnames})

    # ---- LaTeX ----------------------------------------------------------
    # One table per study, with the two schemes as blocks inside it: they share
    # the column structure, and six floats cost more page than three.
    L = []
    present = sorted({r["N"] for r in rows})
    if present:
        L.append(r"\begin{table}[htbp]")
        L.append(r"\centering")
        L.append(r"\small")
        L.append(r"\setlength{\tabcolsep}{4pt}")
        L.append(
            rf"\caption{{Convergence on the {test_name} test, both reconstruction "
            r"schemes. The order is a least-squares fit of $\log E_{\text{shape}}$ "
            r"against $\log N$ over the whole sequence; the residual is the "
            r"largest departure from that fit, in $\log E$. Data: "
            rf"\cite{{figshare2026}}, \protect\path{{{workflow}/results/convergence.csv}}.}}"
        )
        L.append(rf"\label{{tab:{label_prefix}}}")
        L.append(r"\begin{tabular}{r cc cc}")
        L.append(r"\toprule")
        L.append(
            r"& \multicolumn{2}{c}{"
            + FRAME_LABEL[frames[0]]
            + r"} & \multicolumn{2}{c}{"
            + FRAME_LABEL[frames[1]]
            + r"} \\"
        )
        L.append(r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}")
        L.append(
            r"$N$ & $E_{\text{shape}}$ & $E_{\text{bound}}$"
            r" & $E_{\text{shape}}$ & $E_{\text{bound}}$ \\"
        )
        for si, scheme in enumerate(schemes):
            by = {
                frame: {r["N"]: r for r in rows
                        if r["scheme"] == scheme and r["frame"] == frame}
                for frame in frames
            }
            L.append(r"\midrule")
            L.append(rf"\multicolumn{{5}}{{l}}{{\emph{{{scheme}}}}} \\")
            for N in present:
                cells = [str(N)]
                for frame in frames:
                    r = by[frame].get(N)
                    if r is None:
                        cells += ["---", "---"]
                    else:
                        cells += [_fmt(r["E_shape"]), _fmt(r["E_bound"])]
                L.append(" & ".join(cells) + r" \\")
            fo, fr_, _ = fits[(scheme, frames[0])]
            go, gr, _ = fits[(scheme, frames[1])]
            L.append(
                r"\quad fitted order & \multicolumn{2}{c}{"
                + ("---" if fo is None else f"${fo:.2f}$")
                + r"} & \multicolumn{2}{c}{"
                + ("---" if go is None else f"${go:.2f}$")
                + r"} \\"
            )
            L.append(
                r"\quad max residual & \multicolumn{2}{c}{"
                + ("---" if fr_ is None else f"${fr_:.3f}$")
                + r"} & \multicolumn{2}{c}{"
                + ("---" if gr is None else f"${gr:.3f}$")
                + r"} \\"
            )
        L.append(r"\bottomrule")
        L.append(r"\end{tabular}")
        L.append(r"\end{table}")
        L.append("")

    with open(tex_path, "w") as f:
        f.write("\n".join(L))
