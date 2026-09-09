#!/usr/bin/env python3
"""Audit a Cartesian, all-Dirichlet sphere operator exported by PETSc."""
import argparse,csv,json,os,shutil,re
from pathlib import Path
import numpy as np
from scipy.linalg import eigvals
from compare import configure,field,solve
from advanced import mesh_lists
from test_invariants import read_petsc_matrix

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--prepared",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args();args.output=args.output.resolve()
    shutil.copytree(args.prepared,args.output)
    configure(args.output,"transmission")
    os.environ["PETSC_OPTIONS"]="-ksp_view_mat binary:matrix.dat"
    solve(args.output,"acousticHelmholtzFoam",1)
    C=field(args.output/"0/C");V=field(args.output/"0/V",len(C))
    alpha=field(args.output/"0/alpha.water",len(C))
    rho=alpha*1000+(1-alpha)*1.2
    kappa=alpha/(1000*1500**2)+(1-alpha)/(1.2*343**2)
    M=read_petsc_matrix(args.output/"matrix.dat").toarray()[1::2,1::2]
    K=-M/rho[:,None]+np.diag(V*(343/.001)**2*kappa)
    points,faces,_=mesh_lists(args.output)
    text=(args.output/"constant/polyMesh/owner").read_text()
    owner=np.fromstring(re.search(r"\n\s*\d+\s*\n\s*\(\s*\n(.*?)\n\s*\)",text,re.S)[1],sep=" ",dtype=int)
    text=(args.output/"constant/polyMesh/neighbour").read_text()
    ni=int(re.search(r"\n\s*(\d+)\s*\n\s*\(",text)[1])
    neighbour=np.fromstring(re.search(r"\n\s*\d+\s*\n\s*\(\s*\n(.*?)\n\s*\)",text,re.S)[1],sep=" ",dtype=int)
    # Independently integrate the emitted multipoint face operators, retain
    # the original two-point terms elsewhere, and compare with PETSc's matrix.
    expected=np.zeros_like(K)
    replaced=set()
    with (args.output/"postProcessing/acousticInterface/1/operators.tsv").open() as stream:
        for row in csv.DictReader(stream,delimiter="\t"):
            f=int(row["face"]);o=int(row["ownerGlobal"]);n=int(row["neighbourGlobal"])
            k=int(row["columnGlobal"]);w=float(row["weight"])
            expected[o,k]-=w;expected[n,k]+=w;replaced.add(f)
    af=field(args.output/"1/alphaf",ni)
    for f in range(ni):
        if f in replaced: continue
        polygon=points[faces[f]]
        S=.5*np.sum(np.cross(polygon,np.roll(polygon,-1,axis=0)),axis=0)
        o,n=owner[f],neighbour[f]
        g=np.linalg.norm(S)**2/(np.dot(C[n]-C[o],S)*(af[f]*1000+(1-af[f])*1.2))
        expected[o,o]+=g;expected[n,n]+=g;expected[o,n]-=g;expected[n,o]-=g
    KN=K.copy()
    for i in range(ni,len(faces)):
        polygon=points[faces[i]]
        S=.5*np.sum(np.cross(polygon,np.roll(polygon,-1,axis=0)),axis=0)
        d=np.dot(polygon.mean(axis=0)-C[owner[i]],S/np.linalg.norm(S))
        g=np.linalg.norm(S)/(1.2*d)
        KN[owner[i],owner[i]]-=g
        expected[owner[i],owner[i]]+=g
    scale=np.linalg.norm(KN,np.inf)
    ev=eigvals(KN);ed=eigvals(K)
    const=np.linalg.norm(KN@np.ones(len(C)),np.inf)/scale
    zeros=int(np.sum(abs(ev)<1e-10*scale))
    rng=np.random.default_rng(46)
    p=rng.normal(size=len(C))+1j*rng.normal(size=len(C))
    action=np.linalg.norm((K-expected)@p)/(np.linalg.norm(K,np.inf)*np.linalg.norm(p))
    cancellation=np.linalg.norm(np.ones(len(C))@KN,np.inf)/scale
    report={"matrixActionDefect":float(action),"conservationDefect":float(cancellation),"cells":len(C),"constantDefect":float(const),"neumannNullity":zeros,
            "neumannMinRealEigenvalueScaled":float(ev.real.min()/scale),
            "dirichletMinRealEigenvalueScaled":float(ed.real.min()/scale),
            "symmetryDefect":float(np.linalg.norm(KN-KN.T)/np.linalg.norm(KN))}
    report["passed"]=bool(action<=1e-10 and cancellation<=1e-10 and const<=1e-10 and zeros==1 and ev.real.min()>=-1e-10*scale and ed.real.min()>0)
    (args.output/"operator-audit.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    if not report["passed"]: raise SystemExit(1)
if __name__=="__main__":main()
