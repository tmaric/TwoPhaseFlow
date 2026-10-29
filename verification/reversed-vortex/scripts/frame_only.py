#!/usr/bin/env python3
"""The error of the frame alone: alpha = 0 leaves a solid-body revolution.

With the base amplitude switched off the laboratory field is a rigid rotation,
whose exact solution at T is the initial disc. Its error must converge on its
own and stay a small fraction of the framed error; otherwise the frame, not the
cancellation, is what the comparison measures. The runs use the same fixed
time step as the convergence study.
"""

from __future__ import annotations

import argparse
import math
import os

from courant_check import from_log
from frame_identity import series


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--dat", nargs="+", required=True, help="volumeFractionError.dat per N")
    a.add_argument("--ns", required=True, help="comma-separated N, same order as --dat")
    a.add_argument("--end-time", type=float, required=True)
    a.add_argument("--scheme", required=True)
    a.add_argument("--out", required=True)
    args = a.parse_args()

    Ns = [int(n) for n in args.ns.split(",")]
    assert len(Ns) == len(args.dat)
    rows = []
    for N, path in zip(Ns, args.dat):
        s = series(path)
        t = min(s, key=lambda x: abs(x - args.end_time))
        if abs(t - args.end_time) > 1e-6:
            raise SystemExit(f"{path}: no row at T={args.end_time}")
        e_shape, _, _ = s[t]
        e_mass = max(abs(v[1]) for v in s.values())
        _, div, _ = from_log(os.path.join(os.path.dirname(path), "log.advectorVoF"))
        rows.append((N, e_shape, e_mass, div))

    lx = [math.log(r[0]) for r in rows]
    ly = [math.log(r[1]) for r in rows]
    mx, my = sum(lx) / len(lx), sum(ly) / len(ly)
    slope = sum((x - mx) * (y - my) for x, y in zip(lx, ly)) / sum((x - mx) ** 2 for x in lx)

    with open(args.out, "w") as fh:
        fh.write("Frame only (alpha = 0): solid-body revolution of the disc\n")
        fh.write(f"scheme {args.scheme}, end time T = {args.end_time:g}\n\n")
        fh.write(f"{'N':>6} {'E_shape(T)':>14} {'max E_mass':>12} {'max div':>10} {'order':>7}\n")
        for i, (N, e, m, div) in enumerate(rows):
            order = "" if i == 0 else f"{math.log(rows[i-1][1] / e) / math.log(N / rows[i-1][0]):7.2f}"
            ds = f"{div:10.2e}" if div is not None else f"{'-':>10}"
            fh.write(f"{N:6d} {e:14.4e} {m:12.2e} {ds} {order:>7}\n")
        fh.write(f"\nfitted order {-slope:.2f}\n")
    print(open(args.out).read())


if __name__ == "__main__":
    main()
