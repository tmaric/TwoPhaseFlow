#!/usr/bin/env python3
"""Layered interface with PML: original versus geometric interface face areas.

Repeats the published layered configurations of this case with both
interface face-area definitions: both material orderings at NX = 8000, the
five damping values of run_pml_sensitivity.sh and the five aligned meshes of
run_pml_mesh_sensitivity.sh. All other parameters are the caseParams.sh
defaults. Both modes use the same pressure flux (flux legacy) and differ only
in the face area: areaFraction legacy (original donor-based selection) or
plicAverage (geometric averaging). Each configuration is solved serially and
on the requested rank counts; referenceDifference is the relative difference
of the pressure from the serial legacy pressure of the same configuration.
Runs the binaries found in PATH; no build commands.
"""
import argparse, csv, hashlib, json, os, re, shutil, subprocess, time
from pathlib import Path
import numpy as np
from postprocess_interface import analytical_pressure, parse_keyvals

HERE = Path(__file__).resolve().parent
MODES = {"legacy": "legacy", "geometric": "plicAverage"}
CONFIGS = ",".join([f"8000:{s}:forward" for s in (10000, 50000, 100000, 200000, 500000)]
                   + ["8000:500000:reverse"]
                   + [f"{n}:500000:forward" for n in (560, 1120, 2240, 4480, 8960)])

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run(argv, case, log, env=None):
    with (case/log).open("w") as f:
        r = subprocess.run(list(map(str, argv)), cwd=case, env=env, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode:
        raise RuntimeError(f"Exit {r.returncode}: {case/log}")

def field(path, count):
    """Internal field of an ascii volScalarField."""
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
                  Path(path).read_text(), re.S)
    if not m:
        raise ValueError(f"cannot parse {path}")
    v = np.fromstring(m[2], sep=" ")
    if int(m[1]) != count or v.size != count or not np.isfinite(v).all():
        raise ValueError(f"invalid field {path}")
    return v

def prepare(dest, n, sigma, reverse):
    dest.mkdir(parents=True)
    for d in ("0.orig", "0.templates", "constant", "system"):
        shutil.copytree(HERE/d, dest/d)
    shutil.copytree(HERE/"0.orig", dest/"0")
    for f in ("prepareCase", "caseParams.sh"):
        shutil.copy2(HERE/f, dest/f)
    env = dict(os.environ, NX=str(n), SIGMA_MAX=str(sigma))
    if reverse:
        env.update(RHOG="1000", CG="1500", RHOL="1.2", CL="343")
    run(["sh", "./prepareCase"], dest, "log.prepareCase", env)
    for p in (dest/"0").glob("*.in"):
        p.unlink()
    run(["foamDictionary", "system/controlDict", "-entry", "writeFormat", "-set", "ascii"], dest, "log.format")
    run(["foamDictionary", "system/controlDict", "-entry", "writePrecision", "-set", "17"], dest, "log.precision")
    for app in ("blockMesh", "setAlphaField", "setPMLFields"):
        run([app], dest, "log." + app)

def exact(prepared, n):
    """Cell-centre analytical pressure of postprocess_interface.py."""
    p = parse_keyvals(prepared/"studyProperties")
    x = p["X_MIN"] + (np.arange(n) + 0.5)*(p["X_MAX"] - p["X_MIN"])/n
    w = 2*np.pi*p["DRIVE_F"]
    Z1, Z2 = p["RHOG"]*p["CG"], p["RHOL"]*p["CL"]
    return analytical_pressure(x, p["RHOG"], p["CG"], p["CL"], p["PISTON_U"], w/p["CG"], w/p["CL"],
                               p["X_INTERFACE"], p["X_MAX"] - p["PML_L"], p["X_MAX"], -p["SIGMA_MAX"],
                               p["PO"], (Z2 - Z1)/(Z2 + Z1))

def solve(case, exe, mode, ranks):
    """Run the solver; require End, the selected method and residuals <= 1e-10."""
    run(["foamDictionary", "system/fvSchemes", "-entry", "acousticInterface", "-set",
         f"{{ areaFraction {MODES[mode]}; flux legacy; writeDiagnostics false; }}"], case, "log.configure")
    t0 = time.monotonic()
    if ranks > 1:
        run(["foamDictionary", "system/decomposeParDict", "-entry", "numberOfSubdomains", "-set", ranks],
            case, "log.ranks")
        run(["decomposePar"], case, "log.decomposePar")
        run(["mpirun", "-np", ranks, exe, "-parallel"], case, "log.solver")
        run(["reconstructPar", "-latestTime"], case, "log.reconstructPar")
    else:
        run([exe], case, "log.solver")
    seconds = time.monotonic()-t0
    text = (case/"log.solver").read_text()
    if not re.search(r"^End$", text, re.M):
        raise RuntimeError("solver did not finish")
    if f"acousticInterface: areaFraction={MODES[mode]} flux=legacy" not in text:
        raise RuntimeError("selected interface method not confirmed in log")
    res = [float(x) for x in re.findall(r"relativeResidual=([+\deE.-]+)", text)]
    if not res or not np.isfinite(res).all() or max(res) > 1e-10:
        raise RuntimeError(f"residual gate failed: {res}")
    return seconds, res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--configs", default=CONFIGS, help="comma-separated NX:sigmaMax:forward|reverse")
    ap.add_argument("--ranks", default="1,4,8", help="rank counts; the first must be 1")
    args = ap.parse_args()
    ranks = list(map(int, args.ranks.split(",")))
    if ranks[0] != 1:
        raise SystemExit("--ranks must start with 1")
    configs = [(int(n), int(s), o) for n, s, o in (c.split(":") for c in args.configs.split(","))]
    if any(o not in ("forward", "reverse") for _, _, o in configs):
        raise SystemExit("material order must be forward or reverse")
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    exe = Path(shutil.which("acousticHelmholtzFoam")).resolve()
    lib = Path(os.environ["FOAM_USER_LIBBIN"])/"libacousticInterface.so"
    tool = Path(shutil.which("setAlphaField")).resolve()
    hashes = {"acousticHelmholtzFoam": sha(exe), "libacousticInterface": sha(lib), "setAlphaField": sha(tool)}
    def git(*a):
        # None when the case runs from a copy outside the repository
        r = subprocess.run(["git", *a], cwd=HERE, capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None
    (out/"settings.json").write_text(json.dumps(dict(
        vars(args), output=str(out),
        binaries={"solver": str(exe), "setAlphaField": str(tool)}, sha256=hashes,
        scripts={f: sha(HERE/f) for f in ("run_face_area_comparison.py", "postprocess_interface.py",
                                          "prepareCase", "caseParams.sh")},
        revision=git("rev-parse", "HEAD"), caseStatus=git("status", "--porcelain", "."),
        OpenFOAM=os.environ.get("WM_PROJECT_VERSION")), indent=2))
    results = []
    for n, sigma, order in configs:
        tag = f"N{n}-s{sigma}-{order}"
        prepared = out/f"prepared-{tag}"
        prepare(prepared, n, sigma, order == "reverse")
        ref = exact(prepared, n)
        serial = {}
        for mode in MODES:
            for r in ranks:
                case = out/f"{tag}-{mode}-np{r}"
                entry = dict(N=n, sigmaMax=sigma, order=order, mode=mode, ranks=r)
                try:
                    shutil.copytree(prepared, case)
                    seconds, res = solve(case, exe, mode, r)
                    p = field(case/"1/Pre", n) + 1j*field(case/"1/Pim", n)
                    if r == 1:
                        serial[mode] = p
                    diff = lambda q: float(np.linalg.norm(p-q)/np.linalg.norm(q))
                    entry.update(completed=True, seconds=seconds, maxRelativeResidual=max(res),
                                 pressureRelL2=diff(ref), parallelDifference=diff(serial[mode]),
                                 referenceDifference=diff(serial["legacy"]))
                    if entry["referenceDifference"] > 1e-9:
                        raise RuntimeError(f"difference from legacy {entry['referenceDifference']:.3e}")
                except Exception as e:
                    entry.update(completed=False, error=str(e))
                results.append(entry); print(json.dumps(entry), flush=True)
                (out/"results.json").write_text(json.dumps(results, indent=2, allow_nan=False))
    if {"acousticHelmholtzFoam": sha(exe), "libacousticInterface": sha(lib), "setAlphaField": sha(tool)} != hashes:
        raise RuntimeError("Binary changed during the study")
    columns = ["N", "sigmaMax", "order", "mode", "ranks", "pressureRelL2", "referenceDifference",
               "parallelDifference", "seconds", "maxRelativeResidual"]
    done = [r for r in results if r["completed"]]
    with open(out/"face-area-comparison.csv", "w", newline="") as f:
        w = csv.DictWriter(f, columns, extrasaction="ignore")
        w.writeheader(); w.writerows(done)
    if done:
        print(f"runs={len(done)} maxReferenceDifference={max(r['referenceDifference'] for r in done):.3e} "
              f"maxRelativeResidual={max(r['maxRelativeResidual'] for r in done):.3e}")
    failures = len(results) - len(done)
    print("FACE_AREA_COMPARISON_COMPLETED", f"failures={failures}", flush=True)
    if failures: raise SystemExit(1)

if __name__ == "__main__":
    main()
