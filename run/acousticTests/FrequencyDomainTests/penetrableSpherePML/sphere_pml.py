#!/usr/bin/env python3
"""Penetrable sphere with PML truncation in an axisymmetric wedge.

A sphere of radius a = 1 mm on the axis (y) of a one-degree wedge enters only
through setAlphaField volume fractions. An axial point source at y = 2.5a is
an excluded box whose faces carry the exact total pressure. The physical
region r <= 4a, -4a <= y <= 6a is truncated by

  dirichlet     exact total pressure on its outer boundary (no PML),
  pml           Cartesian PML with wedgeAzimuthalStretch (cylindrical PML),
  pmlCartesian  the same PML without the azimuthal stretch.

The host wavenumber is k = 1/a, the cells are uniform with h = a/N, and the
interface method is plicAverage/legacy. Errors are volume-weighted relative
complex-pressure L2 norms in the physical region against sphere_reference.py.
Runs the binaries found in PATH; no build commands.
"""
import argparse, hashlib, json, math, os, re, shutil, subprocess, time
from pathlib import Path
import numpy as np
import sphere_reference

HERE = Path(__file__).resolve().parent
A, D, S = 1e-3, 2.5e-3, 0.25e-3          # sphere radius, source height, hole half-size
K = 1/A                                   # host wavenumber
XP, YMIN, YMAX = 4*A, -4*A, 6*A           # physical region (radius, axial range)
PO = 3
THETA = math.radians(0.5)
# (host rho, c), (drop rho, c); the host takes the gas slot of the solver
MATERIALS = {"waterInAir": ((1.2, 343.0), (1000.0, 1500.0)),
             "oilInWater": ((1000.0, 1500.0), (920.0, 1450.0)),
             "airInWater": ((1000.0, 1500.0), (1.2, 343.0))}

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

def solve(case, exe, ranks):
    """Run the solver; require End, the selected interface method and residuals <= 1e-10."""
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
    if "acousticInterface: areaFraction=plicAverage flux=legacy" not in text:
        raise RuntimeError("selected interface method not confirmed in log")
    res = [float(x) for x in re.findall(r"relativeResidual=([+\deE.-]+)", text)]
    if not res or not np.isfinite(res).all() or max(res) > 1e-10:
        raise RuntimeError(f"residual gate failed: {res}")
    return seconds, res

def exact(points, material):
    (rhog, cg), (rhol, cl) = MATERIALS[material]
    return sphere_reference.pressure(points, A, D, K, K*cg/cl, rhog, rhol)

def physical(C):
    radial = np.hypot(C[:, 0], C[:, 2])
    return (radial <= XP*math.cos(THETA)+1e-12) & (C[:, 1] >= YMIN) & (C[:, 1] <= YMAX)

def prepare(dest, material, truncation, n, sigma_factor, pml_length):
    (rhog, cg), (rhol, cl) = MATERIALS[material]
    freq = K*cg/(2*math.pi)
    pml = truncation != "dirichlet"
    h = A/n
    xmax = XP + (pml_length if pml else 0)
    ylo = YMIN - (pml_length if pml else 0)
    yhi = YMAX + (pml_length if pml else 0)
    nr, ny = round(xmax/h), round((yhi-ylo)/h)
    if abs(nr*h-xmax) > 1e-12 or abs(ny*h-(yhi-ylo)) > 1e-12:
        raise ValueError("mesh not aligned with the physical region")
    for d in ("0", "constant", "system"):
        (dest/d).mkdir(parents=True)
    c, s_ = math.cos(THETA), math.sin(THETA)
    v = [(0, ylo, 0), (xmax*c, ylo, -xmax*s_), (xmax*c, yhi, -xmax*s_), (0, yhi, 0),
         (xmax*c, ylo, xmax*s_), (xmax*c, yhi, xmax*s_)]
    (dest/"system/blockMeshDict").write_text(header("blockMeshDict") + "convertToMeters 1;\nvertices (\n" +
        "\n".join("({:.17g} {:.17g} {:.17g})".format(*p) for p in v) +
        f"\n);\nblocks (hex (0 1 2 3 0 4 5 3) ({nr} {ny} 1) simpleGrading (1 1 1));\nedges ();\n" +
        """boundary (
 axis { type empty; faces ((0 0 3 3)); }
 outer { type patch; faces ((1 2 5 4)); }
 bottom { type patch; faces ((0 1 4 0)); }
 top { type patch; faces ((3 3 5 2)); }
 front { type wedge; faces ((0 3 2 1)); }
 back { type wedge; faces ((0 4 5 3)); }
);
""")
    (dest/"system/topoSetDict").write_text(header("topoSetDict") + f"""actions (
 {{ name keep; type cellSet; action new; source boxToCell; box (-1 {D-S:.17g} -1) ({S:.17g} {D+S:.17g} 1); }}
 {{ name keep; type cellSet; action invert; }}
);
""")
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
acousticInterface { areaFraction plicAverage; flux legacy; writeDiagnostics false; }
""")
    (dest/"system/decomposeParDict").write_text(header("decomposeParDict") + "numberOfSubdomains 8;\nmethod scotch;\n")
    # PMLMin/PMLMax are the outer domain edges; damping acts within L of them.
    sigma = sigma_factor*2*math.pi*freq if pml else 0.0
    (dest/"constant/transportProperties").write_text(header("transportProperties") + f"""
rhol {rhol}; rhog {rhog}; cl {cl}; cg {cg};
kl {1/(rhol*cl**2):.17g}; kg {1/(rhog*cg**2):.17g};
f {freq:.17g};
sigmaMax {sigma:.17g}; po {PO}; PMLType rectangle; L {pml_length:.17g};
wedgeAzimuthalStretch {"true" if truncation == "pml" else "false"};
PMLMin (-100 {ylo:.17g} -100); PMLMax ({xmax:.17g} {yhi:.17g} 100);
""")
    (dest/"system/setAlphaFieldDict").write_text(header("setAlphaFieldDict") +
        f"field alpha.water;\ntype ellipsoid; origin (0 0 0); effectiveRadius {A:.17g}; horizontalLongAxis {A:.17g};\n")
    run(["blockMesh"], dest, "log.blockMesh")
    run(["topoSet"], dest, "log.topoSet")
    run(["subsetMesh", "keep", "-patch", "source", "-overwrite"], dest, "log.subsetMesh")
    # subsetMesh creates the exposed-face patch as 'empty'; it must carry data.
    boundary = "constant/polyMesh/boundary"
    run(["foamDictionary", boundary, "-entry", "entry0/source/type", "-set", "patch"], dest, "log.sourceType")
    run(["foamDictionary", boundary, "-entry", "entry0/source/inGroups", "-remove"], dest, "log.sourceGroups")
    run(["checkMesh"], dest, "log.checkMesh")
    patches = patch_centres(dest)
    if patches["source"][0] != "patch" or not len(patches["source"][1]):
        raise RuntimeError("source patch not usable")
    fields = {"Pre": ("volScalarField", "[1 -1 -2 0 0 0 0]", "0"),
              "Pim": ("volScalarField", "[1 -1 -2 0 0 0 0]", "0"),
              "pa": ("volScalarField", "[1 -1 -2 0 0 0 0]", "0"),
              "alpha.water": ("volScalarField", "[0 0 0 0 0 0 0]", "0"),
              "alphaf": ("surfaceScalarField", "[0 0 0 0 0 0 0]", "0"),
              "rho": ("volScalarField", "[1 -3 0 0 0 0 0]", f"{rhog}"),
              "sigma": ("volTensorField", "[0 0 -1 0 0 0 0]", "(0 0 0 0 0 0 0 0 0)")}
    for nm in ("U", "Ure", "Uim"):
        fields[nm] = ("volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)")
    dirichlet = {"source"} | (set() if pml else {"outer", "top", "bottom"})
    for nm, (kind, dim, zero) in fields.items():
        body = header(nm, kind) + f"dimensions {dim};\ninternalField uniform {zero};\nboundaryField {{\n"
        for patch, (typ, centres) in patches.items():
            if typ in ("empty", "wedge"):
                # surface fields need an explicit value on wedge patches (decomposePar)
                e = f"type {typ};" + (f" value uniform {zero};" if kind.startswith("surface") and typ == "wedge" else "")
            elif nm in ("Pre", "Pim") and patch in dirichlet:
                p = exact(centres, material)
                e = "type fixedValue; value " + values(p.real if nm == "Pre" else p.imag) + ";"
            elif nm in ("alphaf", "sigma"):
                e = f"type calculated; value uniform {zero};"
            else:
                e = "type zeroGradient;"
            body += f"{patch} {{ {e} }}\n"
        (dest/"0"/nm).write_text(body + "}\n")
    run(["setAlphaField"], dest, "log.setAlphaField")
    if pml:
        run(["setPMLFields"], dest, "log.setPMLFields")
    run(["postProcess", "-func", "writeCellCentres", "-time", "0"], dest, "log.centres")
    run(["postProcess", "-func", "writeCellVolumes", "-time", "0"], dest, "log.volumes")
    C = field(dest/"0/C")
    if pml:
        sg = field(dest/"0/sigma", len(C))
        phys = physical(C)
        stretched = "Wedge azimuthal PML stretch" in (dest/"log.setPMLFields").read_text()
        if np.any(sg[phys] != 0) or not np.any(sg[~phys] != 0):
            raise RuntimeError("PML damping not confined to the layer")
        if stretched != (truncation == "pml") or np.any(sg[:, 8] != 0) != stretched:
            raise RuntimeError("azimuthal stretch does not match the truncation")
    return dict(h=h, nr=nr, ny=ny, cells=len(C), frequency=freq, sigmaMax=sigma)

def metrics(case, prepared, material):
    C = field(prepared/"0/C"); V = field(prepared/"0/V", len(C))
    alpha = field(prepared/"0/alpha.water", len(C))
    p = field(case/"1/Pre", len(C)) + 1j*field(case/"1/Pim", len(C))
    ref = exact(C, material)
    phys = physical(C)
    err = abs(p-ref)
    def rel(mask):
        dn = np.sum(V[mask]*abs(ref[mask])**2)
        return float(np.sqrt(np.sum(V[mask]*err[mask]**2)/dn)) if mask.any() and dn > 0 else None
    out = {"physicalRelL2": rel(phys),
           "hostRelL2": rel(phys & (alpha < 1e-8)),
           "dropRelL2": rel(phys & (alpha > 1-1e-8)),
           "interfaceRelL2": rel(phys & (alpha >= 1e-8) & (alpha <= 1-1e-8)),
           "farHostRelL2": rel(phys & (np.linalg.norm(C, axis=1) > 2.5*A)
                               & (np.linalg.norm(C-np.array([0, D, 0]), axis=1) > 2*S)),
           "physicalCells": int(phys.sum())}
    return p[phys], V[phys], out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--materials", default="waterInAir,oilInWater")
    ap.add_argument("--sizes", default="8,16,32,64")
    ap.add_argument("--truncations", default="dirichlet,pml,pmlCartesian")
    ap.add_argument("--ranks", type=int, default=1)
    ap.add_argument("--sigma-factor", type=float, default=12.6, help="sigmaMax/omega")
    ap.add_argument("--pml-thickness", type=float, default=3.25, help="PML thickness in units of a")
    args = ap.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    exe = Path(shutil.which("acousticHelmholtzFoam")).resolve()
    lib = Path(os.environ["FOAM_USER_LIBBIN"])/"libacousticInterface.so"
    tools = {t: Path(shutil.which(t)).resolve() for t in ("setAlphaField", "setPMLFields")}
    hashes = {"acousticHelmholtzFoam": sha(exe), "libacousticInterface": sha(lib),
              **{t: sha(p) for t, p in tools.items()}}
    git = lambda *a: subprocess.check_output(["git", *a], cwd=HERE, text=True).strip()
    pml_length = args.pml_thickness*A
    (out/"settings.json").write_text(json.dumps(dict(
        vars(args), output=str(out), pmlLength=pml_length, sphereRadius=A, sourceHeight=D,
        hostWavenumber=K, po=PO, binaries={t: str(p) for t, p in tools.items()} | {"solver": str(exe)},
        sha256=hashes, scripts={p.name: sha(p) for p in (HERE/"sphere_pml.py", HERE/"sphere_reference.py")},
        revision=git("rev-parse", "HEAD"), caseStatus=git("status", "--porcelain", "."),
        OpenFOAM=os.environ.get("WM_PROJECT_VERSION")), indent=2))
    results, fields = [], {}
    for material in args.materials.split(","):
        (rhog, cg), (rhol, cl) = MATERIALS[material]
        check = sphere_reference.self_check(a=A, d=D, kg=K, cg=cg, cl=cl, rhog=rhog, rhol=rhol)
        (out/material).mkdir()
        (out/material/"reference_check.json").write_text(json.dumps(check, indent=2))
        for n in map(int, args.sizes.split(",")):
            for trunc in args.truncations.split(","):
                prepared = out/material/f"prepared-{trunc}-N{n}"
                case = out/material/f"{trunc}-N{n}"
                entry = dict(material=material, N=n, truncation=trunc, ranks=args.ranks,
                             sigmaOverOmega=args.sigma_factor if trunc != "dirichlet" else 0.0,
                             pmlLength=pml_length if trunc != "dirichlet" else 0.0)
                try:
                    entry.update(prepare(prepared, material, trunc, n, args.sigma_factor, pml_length))
                    shutil.copytree(prepared, case)
                    seconds, residuals = solve(case, exe, args.ranks)
                    p, V, m = metrics(case, prepared, material)
                    fields[(material, trunc, n)] = p
                    entry.update(completed=True, seconds=seconds, maxRelativeResidual=max(residuals), **m)
                    ref = fields.get((material, "dirichlet", n))
                    if trunc != "dirichlet" and ref is not None:
                        entry["minusDirichletRelL2"] = float(np.sqrt(np.sum(V*abs(p-ref)**2)/np.sum(V*abs(ref)**2)))
                except Exception as e:
                    entry.update(completed=False, error=str(e))
                results.append(entry); print(json.dumps(entry), flush=True)
                (out/"results.json").write_text(json.dumps(results, indent=2, allow_nan=False))
    for r in results:
        prev = [q for q in results if q.get("completed") and r.get("completed") and q["material"] == r["material"]
                and q["truncation"] == r["truncation"] and q["N"] == r["N"]//2]
        if prev:
            r["observedOrder"] = math.log2(prev[0]["physicalRelL2"]/r["physicalRelL2"])
    (out/"results.json").write_text(json.dumps(results, indent=2, allow_nan=False))
    if {"acousticHelmholtzFoam": sha(exe), "libacousticInterface": sha(lib),
            **{t: sha(p) for t, p in tools.items()}} != hashes:
        raise RuntimeError("Binary changed during the study")
    failures = sum(not r["completed"] for r in results)
    print("SPHERE_PML_COMPLETED", f"failures={failures}", flush=True)
    if failures: raise SystemExit(1)

if __name__ == "__main__":
    main()
