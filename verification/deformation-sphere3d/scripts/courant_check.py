#!/usr/bin/env python3
"""Confirm the interface Courant number every run actually reached.

The time step is fixed in advance from max||w||_1 (timestep_bound.json), so the
interface Courant number can never exceed max_alpha_co -- unless the bound was
computed for a different case than the one run. That is not hypothetical: the
2D bound was once evaluated with the wrong period and frame. This rule reads
the solver's own per-step report and fails the workflow if any run exceeded
the limit, so the values the manuscript quotes come from the runs themselves.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

LINE = re.compile(r"^Interface Courant Number mean:\s*(\S+)\s+max:\s*(\S+)")


def reached(log):
    peak, steps = 0.0, 0
    with open(log) as fh:
        for line in fh:
            m = LINE.match(line)
            if m:
                steps += 1
                peak = max(peak, float(m.group(2)))
    return peak, steps


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--runs", required=True)
    a.add_argument("--cases", required=True, help="comma-separated case names")
    a.add_argument("--max-alpha-co", type=float, required=True)
    a.add_argument("--tol", type=float, default=1e-3,
                   help="relative tolerance for sampling the bound on a grid")
    a.add_argument("--out", required=True)
    args = a.parse_args()

    limit = args.max_alpha_co * (1.0 + args.tol)
    rows, worst = [], (0.0, None)
    for case in args.cases.split(","):
        peak, steps = reached(os.path.join(args.runs, case, "log.advectorVoF"))
        if steps == 0:
            sys.exit(f"{case}: no interface Courant number in log.advectorVoF")
        rows.append((case, steps, peak))
        if peak > worst[0]:
            worst = (peak, case)

    with open(args.out, "w") as fh:
        fh.write(f"# largest interface Courant number per run; limit {args.max_alpha_co:g}\n")
        fh.write(f"{'case':28s} {'steps':>8s} {'max Co':>10s}\n")
        for case, steps, peak in rows:
            fh.write(f"{case:28s} {steps:8d} {peak:10.4f}\n")
        fh.write(f"\nlargest over all runs: {worst[0]:.4f} ({worst[1]})\n")
    print(open(args.out).read())
    if worst[0] > limit:
        sys.exit(f"interface Courant number {worst[0]:.4f} in {worst[1]} exceeds "
                 f"{args.max_alpha_co:g}: the time-step bound does not match the case")


if __name__ == "__main__":
    main()
