#!/usr/bin/env python3
"""Generate the sphere-translation ladder under Cases/.

Cases/<row>/N<N>/ is a copy of this directory (the base case) with
system/simulationParameter rewritten. Rows (RECONSCHEME, SETALPHAMETHOD,
DAVOFMETHOD):

    davof-quadraticFaces   isoAdvection + davof, the state from initDavofFields
                           (quadratic face fractions)
    davof-detrixheAslam    isoAdvection + davof, the Detrixhe-Aslam state
    plicRDF                isoAdvection + plicRDF, alpha from setAlphaField
    isoAlpha               isoAdvection + isoAlpha, alpha from setAlphaField

N = 32, 64, 128 cells across the unit height (2N along x); the time step is
0.5 h / |U| (Courant 0.5, |U| = 1), fixed, the end time 1.

    python3 genCases.py                       # the full ladder
    python3 genCases.py --N 32 64 --rows davof-quadraticFaces plicRDF
"""
import argparse
import os
import re
import shutil
import stat

HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = {
    "davof-quadraticFaces": {"RECONSCHEME": "davof", "SETALPHAMETHOD": "davof",
                             "DAVOFMETHOD": "quadraticFaces"},
    "davof-detrixheAslam": {"RECONSCHEME": "davof", "SETALPHAMETHOD": "davof",
                            "DAVOFMETHOD": "detrixheAslam"},
    "plicRDF": {"RECONSCHEME": "plicRDF", "SETALPHAMETHOD": "setAlphaField"},
    "isoAlpha": {"RECONSCHEME": "isoAlpha", "SETALPHAMETHOD": "setAlphaField"},
}
NS = [32, 64, 128]
UMAG = 1.0
ENDTIME = 1.0
IGNORE = shutil.ignore_patterns("Cases", "*.py", "*.csv", "*.pdf", "*.png", "*.md",
                                "log.*", "postProcessing", "0", "*.dat",
                                "__pycache__", ".pytest_cache", "*.foam")


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
    ap.add_argument("--N", nargs="+", type=int, default=NS)
    ap.add_argument("--rows", nargs="+", default=list(ROWS))
    ap.add_argument("--out", default=os.path.join(HERE, "Cases"))
    a = ap.parse_args()

    made = []
    for row in a.rows:
        for n in a.N:
            dest = os.path.join(a.out, row, f"N{n}")
            if os.path.isdir(dest):
                shutil.rmtree(dest)
            shutil.copytree(HERE, dest, ignore=IGNORE)
            dt = 0.5*(1.0/n)/UMAG
            vals = {"N": n, "DELTAT": f"{dt:.10g}", "ENDTIME": f"{ENDTIME:g}"}
            vals.update(ROWS[row])
            set_param(os.path.join(dest, "system", "simulationParameter"), vals)
            for s in ("Allrun", "Allclean"):
                p = os.path.join(dest, s)
                os.chmod(p, os.stat(p).st_mode | stat.S_IXUSR | stat.S_IXGRP)
            made.append(dest)
    print(f"[genCases] {len(made)} cases under {a.out}")
    for d in made:
        print("  ", os.path.relpath(d, HERE))


if __name__ == "__main__":
    main()
