#!/usr/bin/env python3
"""Penetrable sphere: original versus geometric interface face areas.

A liquid sphere of radius a = 1 mm (water 1000/1500 in air 1.2/343) is
centred in the cube [-2a, 2a]^3 with N^3 uniform cells and enters only through
setAlphaField volume fractions. An incident plane wave exp(i k x) with
ka = 1 in the gas is scattered by the sphere; the exact total pressure
(regular interior and outgoing exterior spherical waves, orders 0-24) is
prescribed on the cube faces, so no PML is used.

Both modes use the same pressure flux (flux legacy) and differ only in the
interface face area: areaFraction legacy (original donor-based selection) or
plicAverage (geometric averaging). Each mesh is solved serially and on the
requested rank counts; D_np is the volume-weighted relative difference of the
reconstructed parallel pressure from the serial pressure of the same mode.
Runs the binaries found in PATH; no build commands.
"""
import argparse, csv, hashlib, json, os, re, shutil, subprocess, time
from pathlib import Path
import numpy as np
from scipy.special import spherical_jn as jn, spherical_yn as yn, eval_legendre

HERE = Path(__file__).resolve().parent
A, KA = 1e-3, 1.0                         # sphere radius, gas ka
RHOG, CG, RHOL, CL = 1.2, 343.0, 1000.0, 1500.0
FREQ = KA*CG/(2*np.pi*A)
LMAX = 24
MODES = {"legacy": "legacy", "geometric": "plicAverage"}

def header(name, kind="dictionary"):
    return f"FoamFile {{ version 2.0; format ascii; class {kind}; object {name}; }}\n"

def values(x):
    return "nonuniform List<scalar>\n" + str(len(x)) + "\n(\n" + "\n".join(f"{v:.17g}" for v in x) + "\n)"

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run(argv, case, log):
    with (case/log).open("w") as f:
        r = subprocess.run(list(map(str, argv)), cwd=case, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode:
        raise RuntimeError(f"Exit {r.returncode}: {case/log}")

def field(path, count=None):
    """Internal field of an ascii volField; vectors and tensors as rows."""
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s+(\d+)\s*\((.*?)\)\s*;",
                  Path(path).read_text(), re.S)
    if not m:
        raise ValueError(f"cannot parse {path}")
    n, width = int(m[2]), {"scalar": 1, "vector": 3, "tensor": 9}[m[1]]
    v = np.fromstring(m[3].replace("(", " ").replace(")", " "), sep=" ")
    if v.size != n*width or (count is not None and n != count) or not np.isfinite(v).all():
        raise ValueError(f"invalid field {path}")
    return v.reshape(n, width) if width > 1 else v

def patch_centres(case):
    """Type and face centres of every boundary patch."""
    root = case/"constant/polyMesh"
    def body(name):
        t = (root/name).read_text()
        return re.search(r"\n\s*(\d+)\s*\n\s*\(\s*\n(.*?)\n\s*\)\s*\n", t, re.S)[2]
    pts = np.fromstring(body("points").replace("(", " ").replace(")", " "), sep=" ").reshape(-1, 3)
    faces = [np.fromstring(f, sep=" ", dtype=int) for f in re.findall(r"\d+\(([^)]+)\)", body("faces"))]
    out = {}
    for name, raw in re.findall(r"(\w+)\s*\{([^{}]*)\}", (root/"boundary").read_text()):
        if "nFaces" in raw:
            n = int(re.search(r"nFaces\s+(\d+)", raw)[1]); st = int(re.search(r"startFace\s+(\d+)", raw)[1])
            typ = re.search(r"type\s+(\w+)", raw)[1]
            out[name] = (typ, np.array([pts[faces[i]].mean(axis=0) for i in range(st, st+n)]).reshape(-1, 3))
    return out

def exact(points):
    """Total pressure of the incident plane wave and the fluid sphere."""
    kg = KA/A; kl = kg*CG/CL
    r = np.linalg.norm(points, axis=1)
    mu = np.divide(points[:, 0], r, out=np.zeros_like(r), where=r > 0)
    inside = r < A; ext = ~inside
    p = np.zeros(len(r), complex)
    for l in range(LMAX+1):
        inc = (2*l+1)*1j**l
        jg, jg1 = jn(l, kg*A), jn(l, kg*A, True)
        hg, hg1 = jg + 1j*yn(l, kg*A), jg1 + 1j*yn(l, kg*A, True)
        jl, jl1 = jn(l, kl*A), jn(l, kl*A, True)
        M = np.array([[hg, -jl], [kg/RHOG*hg1, -kl/RHOL*jl1]], complex)
        aa, bb = np.linalg.solve(M, np.array([-inc*jg, -inc*kg/RHOG*jg1]))
        P = eval_legendre(l, mu)
        p[inside] += bb*jn(l, kl*r[inside])*P[inside]
        p[ext] += (inc*jn(l, kg*r[ext]) + aa*(jn(l, kg*r[ext]) + 1j*yn(l, kg*r[ext])))*P[ext]
    if not np.isfinite(p).all():
        raise ValueError("non-finite reference pressure")
    return p

def prepare(dest, n):
    for d in ("0", "constant", "system"):
        (dest/d).mkdir(parents=True)
    h = 2*A
    v = [(-h, -h, -h), (h, -h, -h), (h, h, -h), (-h, h, -h),
         (-h, -h, h), (h, -h, h), (h, h, h), (-h, h, h)]
    faces = {"left": "0 4 7 3", "right": "1 2 6 5", "bottom": "0 1 5 4",
             "top": "3 7 6 2", "front": "0 3 2 1", "back": "4 5 6 7"}
    (dest/"system/blockMeshDict").write_text(header("blockMeshDict") + "convertToMeters 1;\nvertices (\n" +
        "\n".join("(%g %g %g)" % p for p in v) +
        f"\n);\nblocks (hex (0 1 2 3 4 5 6 7) ({n} {n} {n}) simpleGrading (1 1 1));\nedges ();\nboundary (\n" +
        "\n".join(f"{p} {{ type patch; faces (({f})); }}" for p, f in faces.items()) + "\n);\n")
    (dest/"system/controlDict").write_text(header("controlDict") + """application acousticHelmholtzFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 1; deltaT 1;
writeControl timeStep; writeInterval 1; writeFormat ascii; writePrecision 17; runTimeModifiable false;
""")
    (dest/"system/fvSolution").write_text(header("fvSolution") + """
solvers { "alpha.*" { reconstructionScheme plicRDF; isoFaceTol 1e-8; surfCellTol 1e-8; } }
SIMPLE { nNonOrthogonalCorrectors 0; }
""")
    (dest/"system/fvSchemes").write_text(header("fvSchemes") + """
ddtSchemes { default steadyState; }
gradSchemes { default leastSquares; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    (dest/"system/decomposeParDict").write_text(header("decomposeParDict") + "numberOfSubdomains 4;\nmethod scotch;\n")
    (dest/"constant/transportProperties").write_text(header("transportProperties") + f"""
rhol {RHOL:g}; rhog {RHOG:g}; cl {CL:g}; cg {CG:g};
kl {1/(RHOL*CL**2):.17g}; kg {1/(RHOG*CG**2):.17g};
f {FREQ:.17g};
sigmaMax 0; po 3; PMLType rectangle; L 1;
PMLMin (-10 -10 -10); PMLMax (10 10 10);
""")
    (dest/"system/setAlphaFieldDict").write_text(header("setAlphaFieldDict") +
        f"field alpha.water;\ntype ellipsoid; origin (0 0 0); effectiveRadius {A:g}; horizontalLongAxis {A:g};\n")
    run(["blockMesh"], dest, "log.blockMesh")
    patches = patch_centres(dest)
    fields = {"Pre": ("volScalarField", "[1 -1 -2 0 0 0 0]", "0"),
              "Pim": ("volScalarField", "[1 -1 -2 0 0 0 0]", "0"),
              "pa": ("volScalarField", "[1 -1 -2 0 0 0 0]", "0"),
              "alpha.water": ("volScalarField", "[0 0 0 0 0 0 0]", "0"),
              "alphaf": ("surfaceScalarField", "[0 0 0 0 0 0 0]", "0"),
              "rho": ("volScalarField", "[1 -3 0 0 0 0 0]", f"{RHOG:g}"),
              "sigma": ("volTensorField", "[0 0 -1 0 0 0 0]", "(0 0 0 0 0 0 0 0 0)")}
    for nm in ("U", "Ure", "Uim"):
        fields[nm] = ("volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)")
    for nm, (kind, dim, zero) in fields.items():
        body = header(nm, kind) + f"dimensions {dim};\ninternalField uniform {zero};\nboundaryField {{\n"
        for patch, (_, centres) in patches.items():
            if nm in ("Pre", "Pim"):
                p = exact(centres)
                e = "type fixedValue; value " + values(p.real if nm == "Pre" else p.imag) + ";"
            elif nm in ("alphaf", "sigma"):
                e = f"type calculated; value uniform {zero};"
            else:
                e = "type zeroGradient;"
            body += f"{patch} {{ {e} }}\n"
        (dest/"0"/nm).write_text(body + "}\n")
    run(["setAlphaField"], dest, "log.setAlphaField")
    run(["postProcess", "-func", "writeCellCentres", "-time", "0"], dest, "log.centres")
    run(["postProcess", "-func", "writeCellVolumes", "-time", "0"], dest, "log.volumes")
    return len(field(dest/"0/V"))

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

def metrics(case, prepared):
    C = field(prepared/"0/C"); V = field(prepared/"0/V", len(C))
    alpha = field(prepared/"0/alpha.water", len(C))
    p = field(case/"1/Pre", len(C)) + 1j*field(case/"1/Pim", len(C))
    ref = exact(C); err = abs(p-ref)
    out = {"pressureRelL2": float(np.sqrt(np.sum(V*err**2)/np.sum(V*abs(ref)**2)))}
    for name, w in (("liquid", alpha), ("gas", 1-alpha),
                    ("interface", ((alpha > 1e-8) & (alpha < 1-1e-8)).astype(float))):
        dn = np.sum(V*w*abs(ref)**2)
        out[name + "RelL2"] = float(np.sqrt(np.sum(V*w*err**2)/dn)) if dn > 0 else None
    return p, V, out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--sizes", default="16,24,32,48")
    ap.add_argument("--modes", default="legacy,geometric")
    ap.add_argument("--ranks", default="1,8", help="rank counts; the first must be 1")
    args = ap.parse_args()
    ranks = list(map(int, args.ranks.split(",")))
    if ranks[0] != 1:
        raise SystemExit("--ranks must start with 1")
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
        vars(args), output=str(out), sphereRadius=A, ka=KA, frequency=FREQ, seriesOrders=LMAX+1,
        binaries={"solver": str(exe), "setAlphaField": str(tool)}, sha256=hashes,
        scripts={"sphere_area_comparison.py": sha(HERE/"sphere_area_comparison.py")},
        revision=git("rev-parse", "HEAD"), caseStatus=git("status", "--porcelain", "."),
        OpenFOAM=os.environ.get("WM_PROJECT_VERSION")), indent=2))
    results = []
    for n in map(int, args.sizes.split(",")):
        prepared = out/f"prepared-N{n}"
        cells = prepare(prepared, n)
        for mode in args.modes.split(","):
            serial = None
            for r in ranks:
                case = out/f"N{n}-{mode}-np{r}"
                entry = dict(N=n, cells=cells, mode=mode, ranks=r)
                try:
                    shutil.copytree(prepared, case)
                    seconds, res = solve(case, exe, mode, r)
                    p, V, m = metrics(case, prepared)
                    if r == 1:
                        serial = p
                    m["parallelDifference"] = float(np.sqrt(np.sum(V*abs(p-serial)**2)/np.sum(V*abs(serial)**2)))
                    if mode == "geometric" and m["parallelDifference"] > 1e-9:
                        raise RuntimeError(f"parallel difference {m['parallelDifference']:.3e}")
                    entry.update(completed=True, seconds=seconds, maxRelativeResidual=max(res), **m)
                except Exception as e:
                    entry.update(completed=False, error=str(e))
                results.append(entry); print(json.dumps(entry), flush=True)
                (out/"results.json").write_text(json.dumps(results, indent=2, allow_nan=False))
    if {"acousticHelmholtzFoam": sha(exe), "libacousticInterface": sha(lib), "setAlphaField": sha(tool)} != hashes:
        raise RuntimeError("Binary changed during the study")
    columns = ["N", "cells", "mode", "ranks", "pressureRelL2", "liquidRelL2", "gasRelL2", "interfaceRelL2",
               "parallelDifference", "seconds", "maxRelativeResidual"]
    with open(out/"area-comparison.csv", "w", newline="") as f:
        w = csv.DictWriter(f, columns, extrasaction="ignore")
        w.writeheader(); w.writerows(r for r in results if r["completed"])
    failures = sum(not r["completed"] for r in results)
    print("SPHERE_AREA_COMPLETED", f"failures={failures}", flush=True)
    if failures: raise SystemExit(1)

if __name__ == "__main__":
    main()
