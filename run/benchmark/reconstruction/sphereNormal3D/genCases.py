#!/usr/bin/env python3
"""Generate the sphere-normal reconstruction ladder under Cases/.

Cases/<scheme>/<setAlphaMethod>/N<N>/ is a copy of this directory (the base case)
with system/simulationParameter rewritten. Standard library only (no casefoam).

    python3 genCases.py                # the full ladder
    python3 genCases.py --N 20 40 --schemes plicRDF gradAlpha --setAlpha cutCellIso
"""
import argparse
import os
import re
import shutil
import stat

HERE = os.path.dirname(os.path.abspath(__file__))
SCHEMES = ["plicRDF", "gradAlpha", "isoAlpha"]
SET_ALPHA = ["cutCellIso", "cutCellImpFunc"]
LADDER = [20, 40, 80, 160]
IGNORE = shutil.ignore_patterns("Cases", "*.py", "*.csv", "*.pdf", "*.png", "*.md",
                                "log.*", "postProcessing", "0", "error.dat",
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
    ap.add_argument("--N", nargs="+", type=int, default=LADDER)
    ap.add_argument("--schemes", nargs="+", default=SCHEMES)
    ap.add_argument("--setAlpha", nargs="+", default=SET_ALPHA)
    ap.add_argument("--out", default=os.path.join(HERE, "Cases"))
    a = ap.parse_args()

    made = []
    for scheme in a.schemes:
        for setA in a.setAlpha:
            for n in a.N:
                dest = os.path.join(a.out, scheme, setA, f"N{n}")
                if os.path.isdir(dest):
                    shutil.rmtree(dest)
                shutil.copytree(HERE, dest, ignore=IGNORE)
                set_param(os.path.join(dest, "system", "simulationParameter"),
                          {"N": n, "RECONSCHEME": scheme, "SETALPHAMETHOD": setA})
                for s in ("Allrun", "Allclean"):
                    p = os.path.join(dest, s)
                    os.chmod(p, os.stat(p).st_mode | stat.S_IXUSR | stat.S_IXGRP)
                made.append(dest)
    print(f"[genCases] {len(made)} cases under {a.out}")
    for d in made:
        print("  ", os.path.relpath(d, HERE))


if __name__ == "__main__":
    main()
