#!/usr/bin/env python3
"""Form the one-step amplification matrix G of acousticWaveFoam by probing.

Usage: amplification.py OUTDIR      (case parameters from the environment,
                                     as for prepareCase; JOBS workers)

The state of every cell is (p, q = pDot, psi = v, w = u) and, for the
segregated update, whose PML pressure equation uses three time levels,
also the previous pressure p^{n-1} (read by OpenFOAM from p_0). Column j
of G is one solver step from e_j + r minus one step from r, with a fixed
random state r: sparse unit states can break down the iterative pressure
solve, shifted ones cannot. G is stored in
OUTDIR/G.npz with the component layout; a step from a random state checks
linearity and the assembly.
"""
import json
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
         "prepare_pml.py")

HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {name};
}}

dimensions      {dims};

internalField   nonuniform List<{kind}> {n}
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
{empty}}}
"""

FIELDS = {
    # name: (class, list kind, dimensions, walls patch entry)
    "p": ("volScalarField", "scalar", "[1 -1 -2 0 0 0 0]",
          "        type            zeroGradient;"),
    "p_0": ("volScalarField", "scalar", "[1 -1 -2 0 0 0 0]",
            "        type            zeroGradient;"),
    "pDot": ("volScalarField", "scalar", "[1 -1 -3 0 0 0 0]",
             "        type            calculated;\n"
             "        value           uniform 0;"),
    "v": ("volScalarField", "scalar", "[1 -1 -1 0 0 0 0]",
          "        type            fixedValue;\n"
          "        value           uniform 0;"),
    "u": ("volVectorField", "vector", "[1 -2 -2 0 0 0 0]",
          "        type            calculated;\n"
          "        value           uniform (0 0 0);"),
}


def write_field(path, name, values, dim):
    cls, kind, dims, wall = FIELDS[name]
    if kind == "scalar":
        lines = "\n".join(repr(float(x)) for x in values)
    else:
        lines = "\n".join("(" + " ".join(repr(float(c)) for c in row) + ")"
                          for row in values)
    empty = ("" if dim == 3 else
             "    frontAndBack\n    {\n        type            empty;\n    }\n")
    (path / name).write_text(HEADER.format(
        cls=cls, name=name, dims=dims, kind=kind, n=len(values),
        values=lines, wall=wall, empty=empty))


def read_field(path, n, vector=False):
    text = path.read_text().split("internalField", 1)[1].split(";", 1)[0]
    if "nonuniform" not in text:
        value = re.search(r"uniform\s+(.*)", text).group(1)
        value = np.array(value.strip("() ").split(), dtype=float)
        return np.tile(value, (n, 1)) if vector else np.full(n, value[0])
    body = text.split("(", 1)[1].rsplit(")", 1)[0]
    data = np.array(body.replace("(", " ").replace(")", " ").split(),
                    dtype=float)
    return data.reshape(n, 3) if vector else data


class Layout:
    """Order of the state vector: component blocks over all cells."""

    def __init__(self, n_cells, dim, three_level):
        self.n, self.dim = n_cells, dim
        self.names = ["p", "q", "psi"] + [f"w{c}" for c in "xyz"[:dim]]
        if three_level:
            self.names.append("pPrev")
        self.size = len(self.names) * n_cells

    def block(self, name):
        i = self.names.index(name)
        return slice(i * self.n, (i + 1) * self.n)

    def to_fields(self, y):
        f = {"p": y[self.block("p")], "pDot": y[self.block("q")],
             "v": y[self.block("psi")]}
        u = np.zeros((self.n, 3))
        for k, c in enumerate("xyz"[:self.dim]):
            u[:, k] = y[self.block(f"w{c}")]
        f["u"] = u
        f["p_0"] = (y[self.block("pPrev")] if "pPrev" in self.names
                    else y[self.block("p")])
        return f

    def from_dir(self, out, p_in):
        y = np.zeros(self.size)
        y[self.block("p")] = read_field(out / "p", self.n)
        y[self.block("q")] = read_field(out / "pDot", self.n)
        y[self.block("psi")] = read_field(out / "v", self.n)
        u = read_field(out / "u", self.n, vector=True)
        for k, c in enumerate("xyz"[:self.dim]):
            y[self.block(f"w{c}")] = u[:, k]
        if "pPrev" in self.names:
            y[self.block("pPrev")] = p_in
        return y


def prepare_base(base):
    if base.exists():
        shutil.rmtree(base)
    base.mkdir(parents=True)
    for name in FILES:
        src = CASE / name
        (shutil.copytree if src.is_dir() else shutil.copy)(src, base / name)
    shutil.copytree(base / "0.orig", base / "0")
    env = dict(os.environ, DOMAIN_PROBES="0")
    # One step: END_TIME = dt, with the defaults of caseParams.sh.
    co, h, c = subprocess.run(
        ["sh", "-c", '. ./caseParams.sh; echo "$COURANT $CELL_SIZE $CG"'],
        cwd=base, env=env, check=True, capture_output=True, text=True
    ).stdout.split()
    env["END_TIME"] = repr(float(co) * float(h) / float(c))
    subprocess.run(["./prepareCase"], cwd=base, env=env, check=True,
                   stdout=subprocess.DEVNULL)
    for app in ("blockMesh", "setPMLFields"):
        with open(base / f"log.{app}", "w") as log:
            subprocess.run([app], cwd=base, check=True, stdout=log,
                           stderr=subprocess.STDOUT)
    control = base / "system" / "controlDict"
    text = control.read_text()
    text = text[:text.index("functions")] + "functions\n{\n}\n"
    text = re.sub(r"writePrecision\s+\d+;", "writePrecision  17;", text)
    control.write_text(text)
    for name in ("p", "p_0", "pDot", "v", "u"):
        (base / "0" / name).unlink(missing_ok=True)
    return json.loads((base / "caseSetup.json").read_text())


def step(work, layout, y):
    """One solver step from state y in directory work; returns the new state."""
    for d in work.iterdir():
        if d.is_dir() and re.fullmatch(r"[0-9.eE+-]+", d.name) and d.name != "0":
            shutil.rmtree(d)
    fields = layout.to_fields(y)
    for name, values in fields.items():
        write_field(work / "0", name, values, layout.dim)
    # Concurrent launches occasionally fail to start the Open MPI singleton
    # (session directory race); such a launch is retried.
    for attempt in range(5):
        with open(work / "log.step", "w") as log:
            result = subprocess.run(["acousticWaveFoam"], cwd=work,
                                    stdout=log, stderr=subprocess.STDOUT)
        if result.returncode == 0:
            break
        if "orte_init failed" not in (work / "log.step").read_text():
            raise RuntimeError(f"{work}: acousticWaveFoam failed")
    else:
        raise RuntimeError(f"{work}: Open MPI did not start in 5 attempts")
    out = [d for d in work.iterdir() if d.is_dir() and d.name != "0"
           and re.fullmatch(r"[0-9.eE+-]+", d.name)]
    if len(out) != 1:
        raise RuntimeError(f"{work}: expected one output time, found {out}")
    return layout.from_dir(out[0], y[layout.block("p")])


def main():
    out = Path(sys.argv[1]).resolve()
    jobs = int(os.environ.get("JOBS", "24"))
    base = out / "base"
    setup = prepare_base(base)
    n_cells = setup["n"] ** setup["dim"]
    three_level = setup["method"] == "thetaNewmark"
    layout = Layout(n_cells, setup["dim"], three_level)

    workers = []
    for k in range(jobs):
        work = out / f"work{k}"
        if work.exists():
            shutil.rmtree(work)
        shutil.copytree(base, work)
        workers.append(work)

    columns = np.zeros((layout.size, layout.size))
    chunks = np.array_split(np.arange(layout.size), jobs)
    rng = np.random.default_rng(0)
    shift = rng.standard_normal(layout.size)
    shifted = step(workers[0], layout, shift)

    def run(k):
        for j in chunks[k]:
            y = shift.copy()
            y[j] += 1.0
            columns[:, j] = step(workers[k], layout, y) - shifted

    with ThreadPoolExecutor(jobs) as pool:
        list(pool.map(run, range(jobs)))

    # Check: a step from a random state against G applied to it.
    rng = np.random.default_rng(1)
    y = rng.standard_normal(layout.size)
    direct = step(workers[0], layout, y)
    check = float(np.linalg.norm(columns @ y - direct) / np.linalg.norm(direct))

    for work in workers:
        shutil.rmtree(work)
    np.savez_compressed(out / "G.npz", G=columns, names=layout.names,
                        nCells=n_cells, superpositionError=check)
    (out / "setup.json").write_text(json.dumps(
        {"setup": setup, "names": layout.names, "nCells": n_cells,
         "superpositionError": check}, indent=2) + "\n")
    print(f"G {layout.size}x{layout.size} ({', '.join(layout.names)}), "
          f"superposition error {check:.2e}")


if __name__ == "__main__":
    main()
