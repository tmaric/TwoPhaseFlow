#!/usr/bin/env python3
"""Verify the selected paper configuration on immutable layered-PML cases."""
import argparse,json,os,re,shutil
from pathlib import Path
import numpy as np
from compare import REPO,configure,field,run,sha,solve,source_snapshot,input_snapshot

def prepare(dest,n,reverse=False,sigma=500000):
    source=REPO/"run/acousticTests/FrequencyDomainTests/layeredInterface1D"
    dest.mkdir(parents=True)
    for name in ["0.orig","0.templates","constant","system"]:
        shutil.copytree(source/name,dest/name)
    shutil.copytree(source/"0.orig",dest/"0")
    for name in ["prepareCase","caseParams.sh"]:shutil.copy2(source/name,dest/name)
    env=os.environ.copy()
    env.update(NX=str(n),CASE_MODE="pml",X_MIN="0",X_MAX=".35",X_INTERFACE=".09",
               PML_XMAX=".35",PML_L=".15",SIGMA_MAX=str(sigma),DRIVE_F="20000",PO="3")
    if reverse:env.update(RHOG="1000",CG="1500",RHOL="1.2",CL="343")
    run(["sh","./prepareCase"],dest,"log.prepare",env)
    for p in (dest/"0").glob("*.in"):p.unlink()
    run(["foamDictionary","system/controlDict","-entry","writeFormat","-set","ascii"],dest,"log.format")
    run(["foamDictionary","system/controlDict","-entry","writePrecision","-set","17"],dest,"log.precision")
    for name in ["blockMesh","setAlphaField","setPMLFields"]:run([name],dest,"log."+name)
    input_snapshot(dest)

def metrics(case,n,reverse,sigma):
    p=field(case/"1/Pre",n)+1j*field(case/"1/Pim",n)
    x=(np.arange(n)+.5)*.35/n
    rho1,c1,rho2,c2=(1000,1500,1.2,343) if reverse else (1.2,343,1000,1500)
    omega=2*np.pi*20000;k1=omega/c1;k2=omega/c2;R=(rho2*c2-rho1*c1)/(rho2*c2+rho1*c1)
    amp=rho1*c1*.01/(1-R*np.exp(2j*k1*.09))
    left=amp*(np.exp(1j*k1*x)+R*np.exp(1j*k1*(.18-x)))
    right=amp*(1+R)*np.exp(1j*(k2*(x-.09)+k1*.09))
    right*=np.exp(-sigma*.15/(4*c2)*(np.maximum(x-.2,0)/.15)**4)
    exact=np.where(x<.09,left,right)
    return p,float(np.linalg.norm(p-exact)/np.linalg.norm(exact))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--baseline",type=Path,required=True)
    ap.add_argument("--ranks",default="1,4")
    ap.add_argument("--remaining-damping",action="store_true",
                    help="Check the three intermediate damping values in the manuscript")
    a=ap.parse_args();a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    a.baseline=a.baseline.resolve()
    exe=Path(shutil.which("acousticHelmholtzFoam")).resolve()
    lib=Path(os.environ["FOAM_USER_LIBBIN"])/"libacousticInterface.so"
    hashes=[sha(exe),sha(lib),sha(a.baseline)]
    source_snapshot(a.output,exe,lib)
    rows=[]
    # Both published 8000-cell checks, the complete aligned mesh sequence,
    # and a weak-damping control. All use the existing published inputs.
    cases=[(8000,False,500000),(8000,True,500000)]+[(n,False,500000) for n in [560,1120,2240,4480,8960]]+[(8000,False,10000)]
    if a.remaining_damping: cases=[(8000,False,s) for s in [50000,100000,200000]]
    for n,reverse,sigma in cases:
        tag=f"N{n}-{'reverse' if reverse else 'forward'}-s{sigma}"
        prepared=a.output/("prepared-"+tag);prepare(prepared,n,reverse,sigma)
        reference=None
        for mode in ["baseline","submission"]:
            for ranks in [int(x) for x in a.ranks.split(",")]:
                case=a.output/f"{tag}-{mode}-np{ranks}";shutil.copytree(prepared,case)
                if mode=="baseline":
                    run(["foamDictionary","system/fvSchemes","-entry","acousticInterface","-remove"],case,"log.baseline")
                # Submission run deliberately inherits the checked-in dictionary.
                seconds,residuals=solve(case,a.baseline if mode=="baseline" else exe,ranks)
                p,error=metrics(case,n,reverse,sigma)
                if mode=="baseline" and ranks==1:reference=p
                if reference is None:raise RuntimeError("Run serial reference first")
                difference=float(np.linalg.norm(p-reference)/np.linalg.norm(reference))
                row=dict(N=n,reverse=reverse,sigma=sigma,mode=mode,ranks=ranks,
                         pressureRelL2=error,baselineDifference=difference,seconds=seconds,residuals=residuals)
                rows.append(row);(a.output/"results.json").write_text(json.dumps(rows,indent=2,allow_nan=False))
                print(json.dumps(row),flush=True)
                if difference>1e-9:raise RuntimeError(f"Submission changed the reference pressure: {case}")
    if [sha(exe),sha(lib),sha(a.baseline)]!=hashes:raise RuntimeError("Binary changed during verification")
    print("SUBMISSION_PML_COMPARISON_PASSED")
if __name__=="__main__":main()
