#!/usr/bin/env python3
"""Generate the reconstruction ladder of this benchmark under Cases/.

Cases/<scheme>/<setAlpha>/N<N>/ is a copy of this directory (the base case) with
system/simulationParameter rewritten. Standard library only (no casefoam). The
same script serves sphereNormal3D and planeNormal3D (the defaults below are
keyed by the directory name).

The setAlpha axis: cutCellIso and cutCellImpFunc (reconstructInterface sets
alpha itself), and davof-<method> (SETALPHAMETHOD davof, DAVOFMETHOD <method>:
initDavofFields writes alpha AND the face fractions alphaf, reconstructInterface
reads both). The davof scheme needs the face fractions, so it runs on the
davof-* states only. The comparators (plicRDF, gradAlpha, isoAlpha) on
davof-quadraticFaces duplicate davof-detrixheAslam (identical cell fractions:
quadraticFaces changes the FACE fractions only) and are skipped unless --all.

    python3 genCases.py                # the default ladder of this benchmark
    python3 genCases.py --N 20 40 --schemes davof plicRDF --setAlpha davof-exactSphere
"""
import argparse
import os
import re
import shutil
import stat

HERE = os.path.dirname(os.path.abspath(__file__))
NAME = os.path.basename(HERE)
SCHEMES = ["davof", "plicRDF", "gradAlpha", "isoAlpha"]
DEFAULTS = {
    "sphereNormal3D": {
        "setAlpha": ["davof-detrixheAslam", "davof-quadraticFaces", "davof-exactSphere",
                     "cutCellIso", "cutCellImpFunc"],
        "N": [20, 40, 80, 160],
    },
    "planeNormal3D": {
        "setAlpha": ["davof-detrixheAslam", "davof-quadraticFaces", "cutCellIso"],
        "N": [20],
    },
}
CFG = DEFAULTS.get(NAME, DEFAULTS["sphereNormal3D"])
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


def params(set_alpha):
    """The SETALPHAMETHOD / DAVOFMETHOD pair of a setAlpha axis value."""
    if set_alpha.startswith("davof-"):
        return {"SETALPHAMETHOD": "davof", "DAVOFMETHOD": set_alpha[len("davof-"):]}
    return {"SETALPHAMETHOD": set_alpha}


def wanted(scheme, set_alpha, all_combos):
    if scheme == "davof":
        return set_alpha.startswith("davof-")
    if set_alpha == "davof-quadraticFaces" and not all_combos:
        return False
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--N", nargs="+", type=int, default=CFG["N"])
    ap.add_argument("--schemes", nargs="+", default=SCHEMES)
    ap.add_argument("--setAlpha", nargs="+", default=CFG["setAlpha"])
    ap.add_argument("--all", action="store_true",
                    help="also the comparators on davof-quadraticFaces")
    ap.add_argument("--out", default=os.path.join(HERE, "Cases"))
    a = ap.parse_args()

    made = []
    for scheme in a.schemes:
        for setA in a.setAlpha:
            if not wanted(scheme, setA, a.all):
                continue
            for n in a.N:
                dest = os.path.join(a.out, scheme, setA, f"N{n}")
                if os.path.isdir(dest):
                    shutil.rmtree(dest)
                shutil.copytree(HERE, dest, ignore=IGNORE)
                vals = {"N": n, "RECONSCHEME": scheme}
                vals.update(params(setA))
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
