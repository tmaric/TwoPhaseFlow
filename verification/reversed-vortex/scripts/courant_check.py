#!/usr/bin/env python3
"""What every run actually reached: Courant number, divergence, volume error.

The time step is fixed in advance from max||w||_1 (timestep_bound.json), so the
interface Courant number can never exceed max_alpha_co -- unless the bound was
computed for a different case than the one run. That is not hypothetical: the
2D bound was once evaluated with the wrong period and frame. This rule reads
the solver's own per-step report and fails the workflow if any run exceeded
the limit.

The same pass collects the per-step discrete divergence max_c|sum_f phi_f|/V_c
that the moving-frame classes print, and the largest volume and boundedness
errors, so every diagnostic the manuscript quotes comes from the runs.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

COURANT = re.compile(r"^Interface Courant Number mean:\s*(\S+)\s+max:\s*(\S+)")
DIVERGENCE = re.compile(r"^movingFrameFlow(?:3D)?: max\|div\(u\)\| = (\S+)")


def from_log(log):
    co, div, steps = 0.0, None, 0
    with open(log) as fh:
        for line in fh:
            m = COURANT.match(line)
            if m:
                steps += 1
                co = max(co, float(m.group(2)))
                continue
            m = DIVERGENCE.match(line)
            if m:
                div = max(div or 0.0, abs(float(m.group(1))))
    return co, div, steps


def from_errors(dat):
    e_mass = e_bound = 0.0
    for line in open(dat):
        r = line.split()
        if not r or r[0].startswith("#"):
            continue
        e_mass = max(e_mass, abs(float(r[2])))
        e_bound = max(e_bound, abs(float(r[3])))
    return e_mass, e_bound


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--runs", required=True)
    a.add_argument("--cases", required=True, help="comma-separated case names")
    a.add_argument("--max-alpha-co", type=float, required=True)
    a.add_argument("--tol", type=float, default=1e-3,
                   help="relative tolerance for sampling the bound on a grid")
    a.add_argument("--out", required=True)
    args = a.parse_args()

    rows = []
    for case in args.cases.split(","):
        d = os.path.join(args.runs, case)
        co, div, steps = from_log(os.path.join(d, "log.advectorVoF"))
        if steps == 0:
            sys.exit(f"{case}: no interface Courant number in log.advectorVoF")
        rows.append((case, steps, co, div) + from_errors(os.path.join(d, "volumeFractionError.dat")))

    def worst(i, keep=lambda r: True):
        sel = [r for r in rows if keep(r) and r[i] is not None]
        return max(sel, key=lambda r: r[i]) if sel else None

    with open(args.out, "w") as fh:
        w = fh.write
        w(f"# per run: steps, largest interface Courant number (limit {args.max_alpha_co:g}),\n"
          "# largest discrete divergence max_c|sum_f phi_f|/V_c (moving-frame runs only),\n"
          "# largest volume error E_mass and boundedness error E_bound over the run\n")
        w(f"{'case':28s} {'steps':>8s} {'max Co':>8s} {'max div':>10s} {'max E_mass':>11s} {'max E_bound':>12s}\n")
        for case, steps, co, div, em, eb in rows:
            ds = f"{div:10.2e}" if div is not None else f"{'-':>10s}"
            w(f"{case:28s} {steps:8d} {co:8.4f} {ds} {em:11.2e} {eb:12.2e}\n")
        w("\nlargest over all runs:\n")
        for label, i in (("interface Courant number", 2), ("discrete divergence", 3),
                         ("E_mass", 4), ("E_bound", 5)):
            r = worst(i)
            if r is not None:
                w(f"  {label:26s} {r[i]:.4g}  ({r[0]}, {r[1]} steps)\n")
    print(open(args.out).read())

    r = worst(2)
    if r[2] > args.max_alpha_co * (1.0 + args.tol):
        sys.exit(f"interface Courant number {r[2]:.4f} in {r[0]} exceeds "
                 f"{args.max_alpha_co:g}: the time-step bound does not match the case")


if __name__ == "__main__":
    main()
