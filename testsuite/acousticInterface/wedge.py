#!/usr/bin/env python3
"""Compare the existing levitator case at fixed meshes and droplet geometries."""
import argparse,json,os,shutil
from pathlib import Path
import numpy as np
from compare import REPO,COMBINATIONS,configure,field,run,sha,solve,source_snapshot,input_snapshot

def prepare(dest,h,aspect):
    source=REPO/"run/acousticTests/FrequencyDomainTests/Andrade2019/WedgeLevitatorDropAlpha"
    dest.mkdir(parents=True)
    for name in ["0.orig","constant","system"]: shutil.copytree(source/name,dest/name)
    for name in ["prepareCase","caseParams.sh","generateCfMeshInputs.py","postMesh"]:
        shutil.copy2(source/name,dest/name)
    for name in ["controlDict","fvSolution","decomposeParDict"]:
        shutil.copy2(Path(__file__).parent/"templates"/(name+".wedge"),dest/"system"/name)
    env=os.environ.copy()
    env.update(RUN_MODE="singlePoint",CFMESH_MAX_CELL_SIZE=str(h),
               DROP_HORIZONTAL_LONG_AXIS=str(.001*aspect),ALPHA_INTERFACE_REFINEMENT="0")
    run(["sh","./prepareCase"],dest,"log.prepare",env)
    for key,value in [("writeFormat","ascii"),("writePrecision","17")]:
        run(["foamDictionary","system/controlDict","-entry",key,"-set",value],dest,"log."+key)
    run(["cartesian2DMesh"],dest,"log.cartesian2DMesh",env)
    run(["extrudeMesh"],dest,"log.extrudeMesh",env)
    run(["sh","./postMesh"],dest,"log.postMesh",env)
    shutil.copytree(dest/"0.orig",dest/"0")
    for p in (dest/"0").glob("*.in"):p.unlink()
    run(["setAlphaField"],dest,"log.setAlphaField",env)
    run(["setPMLFields"],dest,"log.setPMLFields",env)
    run(["postProcess","-func","writeCellCentres","-time","0"],dest,"log.centres")
    run(["postProcess","-func","writeCellVolumes","-time","0"],dest,"log.volumes")
    input_snapshot(dest)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--sizes",default="0.0004,0.0002")
    ap.add_argument("--aspects",default="1,2")
    ap.add_argument("--ranks",default="1,4")
    args=ap.parse_args();args.output=args.output.resolve()
    args.output.mkdir(parents=True,exist_ok=False)
    if not shutil.which("cartesian2DMesh"): raise RuntimeError("cartesian2DMesh is required")
    exe=Path(shutil.which("acousticHelmholtzFoam")).resolve()
    lib=Path(os.environ["FOAM_USER_LIBBIN"])/"libacousticInterface.so"
    hashes=[sha(exe),sha(lib)];source_snapshot(args.output,exe,lib)
    results=[]
    for h in map(float,args.sizes.split(",")):
      for aspect in map(float,args.aspects.split(",")):
        prepared=args.output/f"prepared-h{h}-a{aspect}"
        prepare(prepared,h,aspect)
        C=field(prepared/"0/C");V=field(prepared/"0/V",len(C));reference={}
        for mode in COMBINATIONS:
          for ranks in map(int,args.ranks.split(",")):
            case=args.output/f"h{h}-a{aspect}-{mode}-np{ranks}"
            entry=dict(h=h,aspect=aspect,mode=mode,ranks=ranks,cells=len(V))
            try:
                shutil.copytree(prepared,case);configure(case,mode)
                seconds,residuals=solve(case,exe,ranks)
                p=field(case/"1/Pre",len(V))+1j*field(case/"1/Pim",len(V))
                u=field(case/"1/Ure",len(V))+1j*field(case/"1/Uim",len(V))
                pr=field(case/"1/pr",len(V));field(case/"1/momFlux",len(V))
                force=np.loadtxt(case/"postProcessing/reflectorForce/force.dat",comments="#",ndmin=2)
                if not np.isfinite(force).all(): raise ValueError("Nonfinite reflector force")
                if ranks==1:reference[mode]=p
                difference=float(np.linalg.norm(p-reference[mode])/np.linalg.norm(reference[mode])) if mode in reference else None
                if mode!="legacy" and difference is not None and difference>1e-9:
                    raise RuntimeError(f"Parallel gate failed: {difference}")
                entry.update(completed=True,seconds=seconds,residuals=residuals,
                    pressureRMS=float(np.sqrt(np.sum(V*abs(p)**2)/np.sum(V))),
                    velocityRMS=float(np.sqrt(np.sum(V[:,None]*abs(u)**2)/np.sum(V))),
                    reflectorWedgeForce=float(force[-1,1]),parallelDifference=difference)
            except Exception as e:entry.update(completed=False,error=str(e))
            results.append(entry);print(json.dumps(entry),flush=True)
            (args.output/"results.json").write_text(json.dumps(results,indent=2,allow_nan=False))
    if [sha(exe),sha(lib)]!=hashes: raise RuntimeError("Binary changed during comparison")
    failures=sum(not x.get("completed",False) for x in results)
    print("WEDGE_COMPARISONS_COMPLETED",f"failures={failures}",flush=True)
    if failures: raise SystemExit(1)
if __name__=="__main__":main()
