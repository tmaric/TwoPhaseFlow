#!/usr/bin/env python3
"""Run immutable, dictionary-selected acoustic comparisons. No build commands."""
import argparse, hashlib, json, os, re, shutil, subprocess, tarfile, time
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
COMBINATIONS = {"legacy": ("legacy", "legacy"),
                "geometry": ("plicAverage", "legacy"),
                "transmission": ("plicAverage", "plicTransmission")}

def run(argv, case, log, env=None, expect_failure=False):
    with (case / log).open("w") as output:
        result = subprocess.run(list(map(str, argv)), cwd=case, env=env,
                                stdout=output, stderr=subprocess.STDOUT, timeout=1800)
    if expect_failure:
        if result.returncode == 0:
            raise RuntimeError(f"Expected failure: {case / log}")
    elif result.returncode:
        raise RuntimeError(f"Exit {result.returncode}: {case / log}")
    return result

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def source_snapshot(output,exe,lib):
    paths=["Allwmake","solver/acousticHelmholtzFoam","src/acousticInterface",
           "apps/benchmark/testAcousticInterface","testsuite/acousticInterface"]
    files={}
    for name in paths:
        root=REPO/name
        candidates=[root] if root.is_file() else root.rglob("*")
        for p in candidates:
            if p.is_file() and not p.is_symlink() and "lnInclude" not in p.parts:
                if p.suffix in (".C",".H",".py",".md",".sh",".sbatch") or "templates" in p.parts or "examples" in p.parts or p.name in ("Allwmake","files","options"):
                    files[str(p.relative_to(REPO))]=sha(p)
    with tarfile.open(output/"source-snapshot.tar.gz", "w:gz") as archive:
        for name in sorted(files): archive.add(REPO/name, arcname=name)
    data={"sourceFiles":files,"solver":str(exe),"solverSHA256":sha(exe),
          "librarySHA256":sha(lib),"OpenFOAM":os.environ.get("WM_PROJECT_VERSION"),
          "revision":subprocess.check_output(["git","rev-parse","HEAD"],cwd=REPO,text=True).strip()}
    (output/"source-manifest.json").write_text(json.dumps(data,indent=2))
    (output/"tracked-source.patch").write_bytes(subprocess.check_output(["git","diff","--",*paths],cwd=REPO))

def input_snapshot(prepared):
    files={}
    for dirname in ["0","constant","system"]:
        for p in (prepared/dirname).rglob("*"):
            if p.is_file(): files[str(p.relative_to(prepared))]=sha(p)
    (prepared/"input-manifest.json").write_text(json.dumps(files,indent=2))

def field(path, count=None):
    text = Path(path).read_text()
    match = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if match:
        if count is None: raise ValueError(f"Uniform field needs a cell count: {path}")
        raw=match[1].strip()
        value=np.fromstring(raw.strip("()"), sep=" ") if raw.startswith("(") else float(raw)
        out=np.tile(value, (count, 1)) if np.ndim(value) else np.full(count, value)
    else:
        match=re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s+(\d+)\s*\((.*?)\)\s*;",text,re.S)
        if not match: raise ValueError(f"Cannot parse field: {path}")
        out=np.fromstring(match[3].replace("("," ").replace(")"," "),sep=" ")
        size=int(match[2]); components={"scalar":1,"vector":3,"tensor":9}[match[1]]
        if out.size != size*components: raise ValueError(f"Invalid field length: {path}")
        if components>1: out=out.reshape(size,components)
    if count is not None and len(out)!=count: raise ValueError(f"Wrong cell count: {path}")
    if out.size == 0: raise ValueError(f"No field samples: {path}")
    if not np.isfinite(out).all(): raise ValueError(f"Nonfinite field: {path}")
    return out

def configure(case, mode, diagnostics=True):
    areas,flux=COMBINATIONS[mode]
    d=f"""{{ areaFraction {areas}; flux {flux}; writeDiagnostics {str(diagnostics).lower()};
    plicTransmissionCoeffs {{ maxStencilRings 3; svdRelativeTolerance 1e-12;
    maxConditionNumber 1e10; }} }}"""
    run(["foamDictionary","system/fvSchemes","-entry","acousticInterface","-set",d],
        case,"log.configure")

def prepare_layered(dest,n,offset,reverse=False):
    source=REPO/"run/acousticTests/FrequencyDomainTests/layeredInterface1D"
    dest.mkdir(parents=True)
    for name in ["0.orig","0.templates","constant","system"]:
        shutil.copytree(source/name,dest/name)
    shutil.copytree(source/"0.orig",dest/"0")
    for name in ["prepareCase","caseParams.sh"]:
        shutil.copy2(source/name,dest/name)
    # An isolated prepared case; no Allclean/Allrun or destructive study driver.
    env=os.environ.copy()
    env.update(NX=str(n),CASE_MODE="interface",X_MIN="0",X_MAX="0.14",
               X_INTERFACE=str(0.07+offset*0.14/n),SIGMA_MAX="0",PML_XMAX="0.14")
    if reverse: env.update(RHOG="1000",CG="1500",RHOL="1.2",CL="343")
    run(["sh","./prepareCase"],dest,"log.prepare",env)
    for template in (dest/"0").glob("*.in"): template.unlink()
    # Comparisons start from a neutral snapshot. Submission templates select
    # geometric areas explicitly; each comparison applies its mode afterwards.
    if re.search(r"^\s*acousticInterface\s*\{",(dest/"system/fvSchemes").read_text(),re.M):
        run(["foamDictionary","system/fvSchemes","-entry","acousticInterface","-remove"],dest,"log.neutralInterface")
    run(["foamDictionary","system/controlDict","-entry","writeFormat","-set","ascii"],dest,"log.format")
    run(["foamDictionary","system/controlDict","-entry","writePrecision","-set","17"],dest,"log.precision")
    run(["blockMesh"],dest,"log.blockMesh")
    run(["setAlphaField"],dest,"log.setAlphaField")
    run(["setPMLFields"],dest,"log.setPMLFields")
    return env

def solve(case, executable, ranks):
    start=time.monotonic()
    if ranks>1:
        run(["foamDictionary","system/decomposeParDict","-entry","numberOfSubdomains","-set",ranks],case,"log.ranks")
        run(["decomposePar"],case,"log.decomposePar")
        run(["mpirun","-np",ranks,executable,"-parallel"],case,"log.solver")
        run(["reconstructPar","-latestTime"],case,"log.reconstructPar")
    else: run([executable],case,"log.solver")
    text=(case/"log.solver").read_text()
    if not re.search(r"^End$",text,re.M): raise RuntimeError(f"Incomplete solver: {case}")
    schemes=(case/"system/fvSchemes").read_text()
    area=re.search(r"\bareaFraction\s+(\w+)\s*;",schemes)
    flux=re.search(r"\bflux\s+(\w+)\s*;",schemes)
    if area and flux:
        expected=f"acousticInterface: areaFraction={area[1]} flux={flux[1]}"
        if expected not in text: raise RuntimeError(f"Selected methods not confirmed by solver: {case}")
    residuals=[float(x) for x in re.findall(r"relativeResidual=([+\deE.-]+)",text)]
    # Only the frozen, pre-instrumentation executable may omit this evidence.
    if "acousticInterface:" in text and len(residuals)!=text.count("Acoustic solve:"):
        raise RuntimeError(f"Missing or invalid residual evidence: {case}")
    if "acousticInterface:" in text and not residuals:
        raise RuntimeError(f"No assembled residual reported: {case}")
    if residuals and (not np.isfinite(residuals).all() or max(residuals)>1e-10):
        raise RuntimeError(f"Residual gate failed: {case}: {residuals}")
    changes=[float(x) for x in re.findall(r"Acoustic nonorthogonal change=([+\deE.-]+)",text)]
    if len(changes)>1 and (not np.isfinite(changes).all() or changes[-1]>1e-8):
        raise RuntimeError(f"Nonorthogonal iteration failed: {case}: final change={changes[-1]}")
    return time.monotonic()-start, residuals

def layered_metrics(case,n,offset,reverse):
    p=field(case/"1/Pre",n)+1j*field(case/"1/Pim",n)
    x=(np.arange(n)+0.5)*0.14/n; xi=0.07+offset*0.14/n
    rho1,c1,rho2,c2=(1000,1500,1.2,343) if reverse else (1.2,343,1000,1500)
    # Read frequency actually used by case preparation.
    data=(case/"constant/transportProperties").read_text()
    freq=float(re.search(r"^f\s+(?:\[[^\]]+\]\s+)?([^;]+)",data,re.M)[1])
    k1=2*np.pi*freq/c1; k2=2*np.pi*freq/c2
    R=(rho2*c2-rho1*c1)/(rho2*c2+rho1*c1); T=1+R
    exact=np.where(x<xi,np.exp(1j*k1*(x-xi))+R*np.exp(-1j*k1*(x-xi)),
                   T*np.exp(1j*k2*(x-xi)))
    error=np.abs(p-exact)
    result={"pressureRelL2":float(np.linalg.norm(error)/np.linalg.norm(exact))}
    for name,mask in [("gas",x<xi),("liquid",x>=xi),("interface",abs(x-xi)<2*0.14/n)]:
        result[name+"RelL2"]=float(np.linalg.norm(error[mask])/np.linalg.norm(exact[mask]))
    mask=x<xi-3*0.14/n
    inc,ref=np.linalg.lstsq(np.column_stack([np.exp(1j*k1*(x[mask]-xi)),
                                          np.exp(-1j*k1*(x[mask]-xi))]),p[mask],rcond=None)[0]
    mask=x>xi+3*0.14/n
    transmitted=np.vdot(np.exp(1j*k2*(x[mask]-xi)),p[mask])/mask.sum()
    measuredR,measuredT=ref/inc,transmitted/inc
    result.update(RReal=float(measuredR.real),RImag=float(measuredR.imag),
                  RAbsError=float(abs(measuredR-R)),TAbsError=float(abs(measuredT-T)),
                  RMagnitude=float(abs(measuredR)),RPhase=float(np.angle(measuredR)),
                  TMagnitude=float(abs(measuredT)),TPhase=float(np.angle(measuredT)),
                  RExactMagnitude=float(abs(R)),RExactPhase=float(np.angle(complex(R))),
                  TExactMagnitude=float(abs(T)),TExactPhase=float(np.angle(complex(T))))
    velocity=field(case/"1/Ure",n)+1j*field(case/"1/Uim",n)
    ue=np.where(x<xi,(np.exp(1j*k1*(x-xi))-R*np.exp(-1j*k1*(x-xi)))/(rho1*c1),
                T*np.exp(1j*k2*(x-xi))/(rho2*c2))
    ve=np.linalg.norm(velocity-np.column_stack([ue,np.zeros(n),np.zeros(n)]),axis=1)
    result["velocityRelL2"]=float(np.linalg.norm(ve)/np.linalg.norm(ue))
    alpha=field(case/"0/alpha.water",n)
    for name,mask in [("purePhase0",alpha<1e-8),("purePhase1",alpha>1-1e-8)]:
        result[name+"VelocityRelL2"]=float(np.linalg.norm(ve[mask])/np.linalg.norm(ue[mask]))
    return p,result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--sizes",default="160,240,320,480,640,960")
    ap.add_argument("--offsets",default="0,0.25,0.5")
    ap.add_argument("--ranks",default="1")
    ap.add_argument("--reverse",action="store_true")
    ap.add_argument("--baseline",type=Path)
    args=ap.parse_args()
    args.output=args.output.resolve()
    if args.baseline: args.baseline=args.baseline.resolve()
    args.output.mkdir(parents=True,exist_ok=False)
    exe=Path(shutil.which("acousticHelmholtzFoam")).resolve()
    lib=Path(os.environ["FOAM_USER_LIBBIN"])/"libacousticInterface.so"
    provenance={"revision":subprocess.check_output(["git","rev-parse","HEAD"],cwd=REPO,text=True).strip(),
                "binary":str(exe),"binarySHA256":sha(exe),"librarySHA256":sha(lib),
                "environment":{x:os.environ.get(x) for x in ["WM_PROJECT_VERSION","PETSC_OPTIONS","WM_OPTIONS"]}}
    (args.output/"provenance.json").write_text(json.dumps(provenance,indent=2))
    source_snapshot(args.output,exe,lib)
    results=[]
    for n in map(int,args.sizes.split(",")):
      for offset in map(float,args.offsets.split(",")):
        prepared=args.output/f"prepared-N{n}-o{offset}"
        prepare_layered(prepared,n,offset,args.reverse)
        input_snapshot(prepared)
        modes=list(COMBINATIONS)
        if args.baseline: modes=["baseline","default"]+modes
        reference={}
        for mode in modes:
          for ranks in map(int,args.ranks.split(",")):
            case=args.output/f"N{n}-o{offset}-{mode}-np{ranks}"
            shutil.copytree(prepared,case)
            if mode not in ("baseline","default"): configure(case,mode)
            chosen=args.baseline if mode=="baseline" else exe
            seconds,residuals=solve(case,chosen,ranks)
            p,metrics=layered_metrics(case,n,offset,args.reverse)
            if ranks==1: reference[mode]=p
            if mode in reference:
                metrics["parallelDifference"]=float(np.linalg.norm(p-reference[mode])/max(np.linalg.norm(reference[mode]),1e-300))
                if metrics["parallelDifference"]>1e-9 and mode not in ("baseline","legacy","default"):
                    raise RuntimeError(f"Parallel gate failed: {case}")
            if mode in ("default","legacy") and "baseline" in reference:
                metrics["baselineDifference"]=float(np.linalg.norm(p-reference["baseline"])/np.linalg.norm(reference["baseline"]))
                if metrics["baselineDifference"]>1e-9: raise RuntimeError(f"Legacy gate failed: {case}")
            results.append(dict(N=n,offset=offset,mode=mode,ranks=ranks,seconds=seconds,
                                residuals=residuals,**metrics))
            (args.output/"results.json").write_text(json.dumps(results,indent=2,allow_nan=False))
            print(json.dumps(results[-1]),flush=True)
    if sha(exe)!=provenance["binarySHA256"] or sha(lib)!=provenance["librarySHA256"]:
        raise RuntimeError("Binary/library changed during comparisons")
    print("COMPARISONS_COMPLETED",flush=True)
if __name__=="__main__": main()
