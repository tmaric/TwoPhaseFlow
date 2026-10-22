#!/usr/bin/env python3
"""Amplification matrix of the (p, q) dynamics on the distorted cavity.

Usage: amplification_cavity.py OUTDIR   (case parameters from the environment,
                                         as for prepareCase; JOBS workers)

Without PML the auxiliary fields do not act on the pressure: psi only
integrates p, and w stays zero. The one-step map of (p, q = pDot) is formed
column by column, one solver step from e_j + r minus one step from r for a
fixed random state r, and analysed:
- spectral radius, eigenvalues outside the unit circle;
- for the trapezoidal update the semidiscrete eigenvalues
  lambda = (2/dt)(g - 1)/(g + 1), largest real part, and the corresponding
  amplitude growth per period exp(Re lambda T);
- where the eigenvector of the most unstable mode lives (distance to the
  nearest wall in cells);
- a step from another random state against G (linearity check).
Writes OUTDIR/G.npz and OUTDIR/spectrum.json.
"""
import json
import math
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

CASE = Path(__file__).resolve().parent
FILES = ("0.orig", "constant", "system", "caseParams.sh", "prepareCase",
         "prepare_cavity2d.py")
HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      {name};
}}

dimensions      {dims};

internalField   nonuniform List<scalar> {n}
(
{values}
)
;

boundaryField
{{
    walls
    {{
{wall}
    }}
    frontAndBack
    {{
        type            empty;
    }}
}}
"""
FIELDS = {
    "p": ("[1 -1 -2 0 0 0 0]", "        type            zeroGradient;"),
    "pDot": ("[1 -1 -3 0 0 0 0]", "        type            calculated;\n"
                                  "        value           uniform 0;"),
}


def write_field(path, name, values):
    dims, wall = FIELDS[name]
    (path / name).write_text(HEADER.format(
        name=name, dims=dims, n=len(values), wall=wall,
        values="\n".join(repr(float(x)) for x in values)))


def read_field(path, n):
    text = path.read_text().split("internalField", 1)[1].split(";", 1)[0]
    if "nonuniform" not in text:
        return np.full(n, float(text.split()[-1]))
    return np.array(text.split("(", 1)[1].rsplit(")", 1)[0].split(), dtype=float)


def prepare_base(base):
    if base.exists():
        shutil.rmtree(base)
    base.mkdir(parents=True)
    for name in FILES:
        src = CASE / name
        (shutil.copytree if src.is_dir() else shutil.copy)(src, base / name)
    shutil.copytree(base / "0.orig", base / "0")
    env = dict(os.environ, N_PERIODS="1")
    subprocess.run(["./prepareCase"], cwd=base, env=env, check=True,
                   stdout=subprocess.DEVNULL)
    with open(base / "log.blockMesh", "w") as log:
        subprocess.run(["blockMesh"], cwd=base, check=True, stdout=log,
                       stderr=subprocess.STDOUT)
    subprocess.run(["./prepareCase", "--distort"], cwd=base, env=env,
                   check=True, stdout=subprocess.DEVNULL)
    setup = json.loads((base / "caseSetup.json").read_text())
    control = base / "system/controlDict"
    text = control.read_text()
    text = text[:text.index("functions")] + "functions\n{\n}\n"
    text = re.sub(r"endTime\s+[^;]+;", f"endTime         {setup['deltaT']!r};", text)
    text = re.sub(r"writeInterval\s+[^;]+;", f"writeInterval   {setup['deltaT']!r};", text)
    text = re.sub(r"writePrecision\s+\d+;", "writePrecision  17;", text)
    control.write_text(text)
    for name in ("p", "pDot"):
        (base / "0" / name).unlink(missing_ok=True)
    return setup


def step(work, n, y):
    for d in work.iterdir():
        if d.is_dir() and re.fullmatch(r"[0-9.eE+-]+", d.name) and d.name != "0":
            shutil.rmtree(d)
    write_field(work / "0", "p", y[:n])
    write_field(work / "0", "pDot", y[n:])
    for attempt in range(5):
        with open(work / "log.step", "w") as log:
            result = subprocess.run(["acousticWaveFoam"], cwd=work,
                                    stdout=log, stderr=subprocess.STDOUT)
        if result.returncode == 0:
            break
        if "orte_init failed" not in (work / "log.step").read_text():
            raise RuntimeError(f"{work}: acousticWaveFoam failed")
    out = [d for d in work.iterdir() if d.is_dir() and d.name != "0"
           and re.fullmatch(r"[0-9.eE+-]+", d.name)]
    return np.concatenate((read_field(out[0] / "p", n),
                           read_field(out[0] / "pDot", n)))


def wall_distance(base, setup):
    n = setup["cells"]
    cells = np.arange(n * n)
    i, j = cells % n, cells // n
    return np.minimum(np.minimum(i, n - 1 - i), np.minimum(j, n - 1 - j))


def main():
    out = Path(sys.argv[1]).resolve()
    jobs = int(os.environ.get("JOBS", "24"))
    base = out / "base"
    setup = prepare_base(base)
    n = setup["cells"] ** 2
    size = 2 * n
    workers = []
    for k in range(jobs):
        work = out / f"work{k}"
        if work.exists():
            shutil.rmtree(work)
        shutil.copytree(base, work)
        workers.append(work)

    # The scale of q relative to p is omega; shift and columns in that scale.
    scale = np.concatenate((np.ones(n), np.full(n, setup["omega"])))
    rng = np.random.default_rng(0)
    shift = rng.standard_normal(size) * scale
    shifted = step(workers[0], n, shift)
    g = np.zeros((size, size))
    chunks = np.array_split(np.arange(size), jobs)

    def run(k):
        for j in chunks[k]:
            y = shift.copy()
            y[j] += scale[j]
            g[:, j] = (step(workers[k], n, y) - shifted) / scale[j]

    with ThreadPoolExecutor(jobs) as pool:
        list(pool.map(run, range(jobs)))
    y = rng.standard_normal(size) * scale
    check = float(np.linalg.norm(g @ y - step(workers[0], n, y)) / np.linalg.norm(g @ y))
    for work in workers:
        shutil.rmtree(work)

    # Analysis in the scaled state (p, q/omega).
    gs = g * (1 / scale)[:, None] * scale[None, :]
    values, vectors = np.linalg.eig(gs)
    mod = np.abs(values)
    dt, period = setup["deltaT"], 1 / setup["f"]
    ok = np.abs(values + 1) > 1e-12
    lam = 2 / dt * (values[ok] - 1) / (values[ok] + 1)
    top = int(np.argmax(mod))
    v = np.abs(vectors[:n, top]) ** 2 + np.abs(vectors[n:, top]) ** 2
    v /= v.sum()
    dist = wall_distance(base, setup)
    share = {f"{d}": float(v[dist == d].sum()) for d in range(4)}
    share["4+"] = float(v[dist >= 4].sum())
    unstable = mod > 1 + 1e-10
    result = {
        "setup": setup, "superpositionError": check,
        "spectralRadius": float(mod.max()), "nUnstable": int(unstable.sum()),
        "topEigenvalue": [float(values[top].real), float(values[top].imag)],
        "amplitudeGrowthPerPeriodFromG": float(mod.max() ** setup["stepsPerPeriod"] - 1),
        "semidiscreteMaxRe": float(lam.real.max()),
        "semidiscreteMaxReRelative": float(lam.real.max() / np.abs(lam).max()),
        "amplitudeGrowthPerPeriodSemidiscrete": float(math.exp(max(lam.real.max(), 0) * period) - 1),
        "semidiscretePositiveRe": int(np.sum(lam.real > 1e-9 * np.abs(lam).max())),
        "topModeShareByWallDistance": share,
        "topModeFrequencyOverMode": float(abs(np.angle(values[top])) / dt / setup["omega"]),
    }
    np.savez_compressed(out / "G.npz", G=g)
    (out / "spectrum.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"{out.name}: rho - 1 = {result['spectralRadius'] - 1:+.3e} "
          f"({result['nUnstable']} outside), growth/period from G "
          f"{result['amplitudeGrowthPerPeriodFromG']:+.3e}; semidiscrete max Re "
          f"{result['semidiscreteMaxRe']:+.3e} 1/s ({result['semidiscretePositiveRe']} "
          f"with Re > 0), growth/period {result['amplitudeGrowthPerPeriodSemidiscrete']:+.3e}; "
          f"top mode by wall distance {share}; superposition {check:.1e}")


if __name__ == "__main__":
    main()
