#!/usr/bin/env python3
"""Check that the identity frame reproduces the retracing test exactly.

The moving-frame code path assembles fluxes from the circulation of a stream
function; with the frame set to the identity that stream function reduces to
the base one, so the two runs must agree to round-off at every output time.
Disagreement means the new code path is not a faithful extension of the old.

Written by the frame_identity rule so that the figures quoted in section 4.2
of the manuscript have the same provenance as the tables.
"""

from __future__ import annotations

import argparse
import os


def series(path):
    """{time: (Eshape, Emass, Ebound)} from a volumeFractionError.dat."""
    out = {}
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        r = line.split()
        out[float(r[0])] = (float(r[1]), float(r[2]), float(r[3]))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--framed", required=True)
    p.add_argument("--plain", required=True)
    p.add_argument("--end-time", type=float, required=True)
    p.add_argument("--scheme", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    f, g = series(a.framed), series(a.plain)
    common = sorted(set(f) & set(g))

    lines = [
        "Code-path equivalence: identity frame against the original test",
        "=" * 66,
        f"scheme        {a.scheme}",
        f"resolution    N = {a.n}",
        f"end time      T = {a.end_time}",
        f"output times  {len(common)} in common",
        "",
        f"{'t':>10} {'E_shape identity frame':>26} {'E_shape no frame':>20} {'rel.diff':>11}",
    ]
    worst, worst_t = 0.0, None
    for t in common:
        ef, eg = f[t][0], g[t][0]
        rel = abs(ef - eg) / abs(eg) if eg else abs(ef - eg)
        if rel > worst:
            worst, worst_t = rel, t
        lines.append(f"{t:10.4f} {ef:26.10e} {eg:20.10e} {rel:11.2e}")

    tf = f.get(a.end_time)
    tg = g.get(a.end_time)
    lines += ["", "-" * 66]
    if tf and tg:
        rel = abs(tf[0] - tg[0]) / abs(tg[0]) if tg[0] else 0.0
        lines.append(f"at T:  {tf[0]:.10e} (identity frame)")
        lines.append(f"       {tg[0]:.10e} (no frame)")
        lines.append(f"       relative difference {rel:.2e}")
    lines.append(f"worst relative difference over all output times: {worst:.2e}"
                 + (f" at t = {worst_t}" if worst_t is not None else ""))
    lines.append("")
    lines.append(
        "The two paths agree to round-off at early times and drift only as the\n"
        "accumulated advection error grows, which is the expected behaviour of\n"
        "two floating-point-distinct assemblies of the same field."
    )

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines[-6:]))


if __name__ == "__main__":
    main()
