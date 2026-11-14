#!/usr/bin/env python3
"""Regenerate the convergence results of the three studies from the run data.

Needs neither OpenFOAM nor SLURM: it reads only

    verification/<study>/config.yaml
    verification/<study>/runs/<case>/volumeFractionError.dat

which the data archive contains, and calls the same functions the Snakemake
workflows call (scripts/aggregate.py, scripts/plot_convergence.py). It writes

    convergence.csv, convergence_table.tex, convergence_shape_error.{png,pdf}

per study into a separate directory, so the result can be compared with the
archived verification/<study>/results/:

    python3 verification/postprocess.py                  # -> results-regenerated/
    python3 verification/postprocess.py --out DIR
    diff verification/reversed-vortex/results/convergence_table.tex \\
         results-regenerated/reversed-vortex/convergence_table.tex

The interface snapshots need the volume-fraction fields, which the archive
leaves out for size; rerun the study (./run.sh) to regenerate them.
"""
import argparse
import importlib.util
import os
import subprocess
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))

# study directory -> title of the convergence figure (as in each Snakefile)
STUDIES = {
    "reversed-vortex": "Reversed vortex, 2D",
    "deformation-sphere3d": "3D shear flow",
    "leveque3d": "LeVeque deformation, 3D",
}


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--out", default=os.path.join(os.getcwd(), "results-regenerated"))
    p.add_argument("studies", nargs="*", default=list(STUDIES))
    a = p.parse_args()

    for study in a.studies:
        wf = os.path.join(HERE, study)
        with open(os.path.join(wf, "config.yaml")) as f:
            cfg = yaml.safe_load(f)
        out = os.path.join(os.path.abspath(a.out), study)
        os.makedirs(out, exist_ok=True)

        aggregate = _load(os.path.join(wf, "scripts", "aggregate.py"), f"aggregate_{study}")
        aggregate.write(
            runs=os.path.join(wf, "runs"),
            schemes=list(cfg["schemes"]),
            frames=list(cfg["frames"]),
            resolutions=[int(n) for n in cfg["resolutions"]],
            end_time=cfg["end_time"],
            csv_path=os.path.join(out, "convergence.csv"),
            tex_path=os.path.join(out, "convergence_table.tex"),
        )
        subprocess.run(
            [sys.executable, os.path.join(wf, "scripts", "plot_convergence.py"),
             "--csv", os.path.join(out, "convergence.csv"),
             "--out", os.path.join(out, "convergence_shape_error.png"),
             "--title", STUDIES[study]],
            check=True,
        )
        print(f"{study}: wrote {out}")


if __name__ == "__main__":
    main()
