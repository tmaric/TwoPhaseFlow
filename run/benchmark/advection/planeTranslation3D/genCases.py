#!/usr/bin/env python3
"""Generate the plane-translation gate cases under Cases/.

Cases/<velocity>/<davofMethod>/N<N>/ is a copy of this directory (the base
case) with system/simulationParameter rewritten: the velocity (normal: U = n,
the plane's unit normal; oblique: U = (1 0 0)), the DAVOF state construction
(detrixheAslam, quadraticFaces: both exact for a plane) and N. The time step
is 0.5 h / |U| (|U| = 1 for both velocities) and the end time ten steps.

    python3 genCases.py
    python3 genCases.py --N 20 40 --velocity normal --davofMethod quadraticFaces
"""
import argparse
import os
import re
import shutil
import stat

HERE = os.path.dirname(os.path.abspath(__file__))
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
    ap.add_argument("--N", nargs="+", type=int, default=[20])
    ap.add_argument("--velocity", nargs="+", default=list(VELOCITIES))
    ap.add_argument("--davofMethod", nargs="+", default=METHODS)
    ap.add_argument("--out", default=os.path.join(HERE, "Cases"))
    a = ap.parse_args()

    base = read_param(os.path.join(HERE, "system", "simulationParameter"))
    L = float(base["L"])
    made = []
    for vel in a.velocity:
        ux, uy, uz = VELOCITIES[vel]
        umag = (ux*ux + uy*uy + uz*uz)**0.5
        for method in a.davofMethod:
            for n in a.N:
                dest = os.path.join(a.out, vel, method, f"N{n}")
                if os.path.isdir(dest):
                    shutil.rmtree(dest)
                shutil.copytree(HERE, dest, ignore=IGNORE)
                dt = 0.5*(L/n)/umag
                set_param(os.path.join(dest, "system", "simulationParameter"),
                          {"N": n, "VELOCITY": vel, "UX": ux, "UY": uy, "UZ": uz,
                           "DELTAT": f"{dt:.10g}", "ENDTIME": f"{NSTEPS*dt:.10g}",
                           "DAVOFMETHOD": method})
                for s in ("Allrun", "Allclean"):
                    p = os.path.join(dest, s)
                    os.chmod(p, os.stat(p).st_mode | stat.S_IXUSR | stat.S_IXGRP)
                made.append(dest)
    print(f"[genCases] {len(made)} cases under {a.out}")
    for d in made:
        print("  ", os.path.relpath(d, HERE))


if __name__ == "__main__":
    main()
