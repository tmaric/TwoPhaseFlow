#!/usr/bin/env python3
"""Generate the plane-translation gate cases under Cases/.

Cases/<velocity>/<row>/N<N>/ is a copy of this directory (the base case)
with system/simulationParameter rewritten: the velocity (normal: U = n, the
plane's unit normal; oblique: U = (1 0 0)), the row (../advectionRows.py:
the davof scheme on the initDavofFields state with its normal model
normalModel/applyTo/sliverTol/maxAngle, the ring rows, and the comparators
plicRDF and isoAlpha on the setAlphaField state) and N. The time step is
Co h / |U| (--Co, default 0.2; |U| = 1 for both velocities; Scheufler and
Roenby 2019 demonstrate second order for CFL 0.2 and below) and the end
time ten steps; Co is written into the case (CO) and the results CSV.

Default: the cluster list, velocity normal (the arm that grew), N = 20 and
40, the state detrixheAslam (quadraticFaces is identical for a plane), rows
davof, Youngs 0.01 / 0.1 / 1, exact, plicRDF, the ring rows ringArea,
ringArea-check, ringFit, ringFit-check, and the comparators.

    python3 genCases.py                                   # the cluster list
    python3 genCases.py --N 20 --velocity normal oblique --rows davof-detrixheAslam
    python3 genCases.py --Co 0.5                          # the time step of the 2026-10-13 cluster runs
    python3 genCases.py --list                            # print the case list only
"""
import argparse
import os
import re
import shutil
import stat
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import advectionRows  # noqa: E402

VELOCITIES = {"normal": (0.3, 0.5, 0.8124), "oblique": (1.0, 0.0, 0.0)}
METHODS = ["detrixheAslam", "quadraticFaces"]
NSTEPS = 10
IGNORE = shutil.ignore_patterns("Cases", "*.py", "*.csv", "*.pdf", "*.png", "*.md",
                                "log.*", "postProcessing", "0", "*.dat",
                                "__pycache__", ".pytest_cache", "*.foam")


def read_param(path):
    vals = {}
    for line in open(path):
        m = re.match(r"^\s*([A-Za-z_]\w*)\s+([^;]+);", line)
        if m:
            vals[m.group(1)] = m.group(2).strip()
    return vals


def set_param(path, values):
    """Rewrite `KEY value;` lines of system/simulationParameter."""
    with open(path) as fh:
        text = fh.read()
    for key, val in values.items():
        text, n = re.subn(rf"^{key}\s+[^;]*;", f"{key} {val};", text, flags=re.M)
        if n != 1:
            raise RuntimeError(f"{path}: expected one line for {key}, found {n}")
    with open(path, "w") as fh:
        fh.write(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--N", nargs="+", type=int, default=[20, 40])
    ap.add_argument("--velocity", nargs="+", default=["normal"], choices=list(VELOCITIES))
    ap.add_argument("--davofMethod", nargs="+", default=["detrixheAslam"], choices=METHODS)
    ap.add_argument("--youngsTol", nargs="*", type=float, default=advectionRows.DEFAULT_YOUNGS_TOLS,
                    help="sliverTol values of the Youngs rows")
    ap.add_argument("--sliverPlicRDFTol", nargs="*", type=float, default=[],
                    help="sliverTol values of the sliver-plicRDF rows (none by default)")
    ap.add_argument("--rows", nargs="+", default=None,
                    help="row names to generate (default: every row of the list)")
    ap.add_argument("--Co", type=float, default=0.2,
                    help="Courant number of the fixed time step Co h / |U| (default 0.2)")
    ap.add_argument("--noComparators", action="store_true")
    ap.add_argument("--list", action="store_true", help="print the case list, generate nothing")
    ap.add_argument("--out", default=os.path.join(HERE, "Cases"))
    a = ap.parse_args()

    rows = advectionRows.all_rows(a.davofMethod, a.youngsTol, a.sliverPlicRDFTol,
                                  comparators=not a.noComparators)
    if a.rows:
        unknown = [r for r in a.rows if r not in rows]
        if unknown:
            sys.exit(f"[genCases] unknown rows {unknown}; known: {list(rows)}")
        rows = {r: rows[r] for r in a.rows}

    base = read_param(os.path.join(HERE, "system", "simulationParameter"))
    L = float(base["L"])
    made = []
    for vel in a.velocity:
        ux, uy, uz = VELOCITIES[vel]
        umag = (ux*ux + uy*uy + uz*uz)**0.5
        for row, vals in rows.items():
            for n in a.N:
                dest = os.path.join(a.out, vel, row, f"N{n}")
                if a.list:
                    print(os.path.relpath(dest, HERE))
                    continue
                if os.path.isdir(dest):
                    shutil.rmtree(dest)
                shutil.copytree(HERE, dest, ignore=IGNORE)
                dt = a.Co*(L/n)/umag
                p = {"N": n, "VELOCITY": vel, "UX": ux, "UY": uy, "UZ": uz, "CO": f"{a.Co:g}",
                     "DELTAT": f"{dt:.10g}", "ENDTIME": f"{NSTEPS*dt:.10g}"}
                p.update(advectionRows.complete(vals))
                set_param(os.path.join(dest, "system", "simulationParameter"), p)
                for s in ("Allrun", "Allclean"):
                    sp = os.path.join(dest, s)
                    os.chmod(sp, os.stat(sp).st_mode | stat.S_IXUSR | stat.S_IXGRP)
                made.append(dest)
    if a.list:
        return
    print(f"[genCases] {len(made)} cases under {a.out}")
    for d in made:
        print("  ", os.path.relpath(d, HERE))


if __name__ == "__main__":
    main()
