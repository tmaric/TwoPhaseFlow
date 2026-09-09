#!/usr/bin/env python3
"""Analytical waves, oblique mesh cuts, and penetrable-sphere comparisons."""
import argparse,json,math,os,re,shutil,subprocess,time
from pathlib import Path
import numpy as np
from scipy.special import spherical_jn as jn, spherical_yn as yn, eval_legendre
from compare import REPO,COMBINATIONS,configure,field,run,sha,solve,source_snapshot,input_snapshot

def header(name,kind="dictionary",location=None):
    return f'FoamFile {{ version 2.0; format ascii; class {kind}; object {name}; }}\n'
def values(x):
    x=np.asarray(x)
    return "nonuniform List<scalar>\n"+str(len(x))+"\n(\n"+"\n".join(f"{v:.17g}" for v in x)+"\n)"
def mesh_lists(case):
    root=case/"constant/polyMesh"
    def body(name):
        s=(root/name).read_text()
        return re.search(r"\n\s*(\d+)\s*\n\s*\(\s*\n(.*?)\n\s*\)",s,re.S)[2]
    points=np.fromstring(body("points").replace("("," ").replace(")"," "),sep=" ").reshape(-1,3)
    faces=[np.fromstring(s,sep=" ",dtype=int) for s in re.findall(r"\d+\(([^)]+)\)",body("faces"))]
    patches={}
    for name,raw in re.findall(r"(\w+)\s*\{([^{}]*)\}",(root/"boundary").read_text()):
        if "nFaces" in raw:
            count=int(re.search(r"nFaces\s+(\d+)",raw)[1]);start=int(re.search(r"startFace\s+(\d+)",raw)[1])
            typ=re.search(r"type\s+(\w+)",raw)[1]
            patches[name]=(typ,np.array([points[faces[i]].mean(axis=0) for i in range(start,start+count)]))
    return points,faces,patches

def plane_pressure(points,omega=2*np.pi*20000,xi=0.0001,homogeneous=False):
    x=points[:,0]-xi
    if homogeneous: return np.exp(1j*omega/343*x)
    R=(1000*1500-1.2*343)/(1000*1500+1.2*343)
    return np.where(x<0,np.exp(1j*omega/343*x)+R*np.exp(-1j*omega/343*x),
                    (1+R)*np.exp(1j*omega/1500*x))

def sphere_pressure(points,a=0.001,ka=1.0,aspect=1.0):
    """Independent regular/outgoing spherical-wave solution; liquid interior."""
    kg=ka/a;kl=kg*343/1500
    r=np.linalg.norm(points,axis=1)
    mu=np.divide(points[:,0],r,out=np.zeros_like(r),where=r>0)
    inside=r<a; p=np.zeros(len(r),complex)
    for ell in range(25):
        inc=(2*ell+1)*1j**ell
        jg=jn(ell,kg*a);jg1=jn(ell,kg*a,True)
        hg=jg+1j*yn(ell,kg*a);hg1=jg1+1j*yn(ell,kg*a,True)
        jl=jn(ell,kl*a);jl1=jn(ell,kl*a,True)
        mat=np.array([[hg,-jl],[kg/1.2*hg1,-kl/1000*jl1]],complex)
        aa,bb=np.linalg.solve(mat,np.array([-inc*jg,-inc*kg/1.2*jg1]))
        leg=eval_legendre(ell,mu)
        p[inside]+=bb*jn(ell,kl*r[inside])*leg[inside]
        ext=~inside
        p[ext]+=(inc*jn(ell,kg*r[ext])+aa*(jn(ell,kg*r[ext])+1j*yn(ell,kg*r[ext])))*leg[ext]
    if not np.isfinite(p).all(): raise ValueError("Invalid analytical sphere pressure")
    return p

def prepare(case,scenario,n,angle=0,aspect=1,mesh_style="perturbed"):
    case.mkdir(parents=True)
    for d in ["0","constant","system"]: (case/d).mkdir()
    sphere=scenario in ("sphere","ellipsoid","sphereWedge")
    xmin,xmax=(-.002,.002) if sphere else (-.07,.07)
    yh=zh=.002 if sphere else .001
    ny=n if sphere else 8
    nz=n if sphere else 1
    xyz=[(xmin,-yh,-zh),(xmax,-yh,-zh),(xmax,yh,-zh),(xmin,yh,-zh),
         (xmin,-yh,zh),(xmax,-yh,zh),(xmax,yh,zh),(xmin,yh,zh)]
    faceverts={"left":"0 4 7 3","right":"1 2 6 5","bottom":"0 1 5 4",
               "top":"3 7 6 2","front":"0 3 2 1","back":"4 5 6 7"}
    patchtypes={p:("patch" if sphere or p in ("left","right") else
                  "empty" if p in ("front","back") else "symmetryPlane") for p in faceverts}
    mesh=header("blockMeshDict")+ "convertToMeters 1;\nvertices (\n"
    mesh+="\n".join("(%g %g %g)"%v for v in xyz)+f"\n);\nblocks (hex (0 1 2 3 4 5 6 7) ({n} {ny} {nz}) simpleGrading (1 1 1));\nedges ();\nboundary (\n"
    mesh+="\n".join(f"{p} {{ type {patchtypes[p]}; faces (({v})); }}" for p,v in faceverts.items())+"\n);\n"
    if scenario=="sphereWedge":
        radius=.002;half=.002;theta=np.deg2rad(.5)
        v=[(0,-half,0),(radius*np.cos(theta),-half,-radius*np.sin(theta)),
           (radius*np.cos(theta),half,-radius*np.sin(theta)),(0,half,0),
           (radius*np.cos(theta),-half,radius*np.sin(theta)),
           (radius*np.cos(theta),half,radius*np.sin(theta))]
        mesh=header("blockMeshDict")+"convertToMeters 1;\nvertices (\n"
        mesh+="\n".join("({:.17g} {:.17g} {:.17g})".format(*x) for x in v)
        mesh+=f"\n);\nblocks (hex (0 1 2 3 0 4 5 3) ({n} {2*n} 1) simpleGrading (1 1 1));\n"
        mesh+="""edges (); boundary (
 axis { type empty; faces ((0 0 3 3)); }
 outer { type patch; faces ((1 2 5 4)); }
 bottom { type patch; faces ((0 1 4 0)); }
 top { type patch; faces ((3 3 5 2)); }
 front { type wedge; faces ((0 3 2 1)); }
 back { type wedge; faces ((0 4 5 3)); }
);\n"""
    (case/"system/blockMeshDict").write_text(mesh)
    (case/"system/controlDict").write_text(header("controlDict")+"""application acousticHelmholtzFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 1; deltaT 1;
writeControl timeStep; writeInterval 1; writeFormat ascii; writePrecision 17; runTimeModifiable false;
""")
    (case/"system/fvSolution").write_text(header("fvSolution")+f"""
solvers {{ "alpha.*" {{ reconstructionScheme plicRDF; isoFaceTol 1e-8; surfCellTol 1e-8; }} }}
SIMPLE {{ nNonOrthogonalCorrectors {100 if angle else 0}; }}
""")
    (case/"system/fvSchemes").write_text(header("fvSchemes")+"""
ddtSchemes { default steadyState; }
gradSchemes { default leastSquares; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    (case/"system/decomposeParDict").write_text(header("decomposeParDict")+"numberOfSubdomains 4;\nmethod scotch;\n")
    frequency=343/(2*np.pi*.001) if sphere else 2000
    (case/"constant/transportProperties").write_text(header("transportProperties")+f"""
rhol 1000; rhog 1.2; cl 1500; cg 343;
kl {1/(1000*1500**2):.17g}; kg {1/(1.2*343**2):.17g};
f {frequency:.17g};
sigmaMax 0; po 3; PMLType rectangle; L 1;
PMLMin (-10 -10 -10); PMLMax (10 10 10);
""")
    run(["blockMesh"],case,"log.blockMesh")
    if angle:
        root=case/"constant/polyMesh/points";text=root.read_text()
        match=re.search(r"(\n\s*\d+\s*\n\s*\(\s*\n)(.*?)(\n\s*\))",text,re.S)
        points=np.fromstring(match[2].replace("("," ").replace(")"," "),sep=" ").reshape(-1,3)
        shape=np.sin(np.pi*(points[:,0]-xmin)/(xmax-xmin)) if mesh_style=="perturbed" else 1
        points[:,0]+=np.tan(np.deg2rad(angle))*points[:,1]*shape
        replacement="\n".join("({:.17g} {:.17g} {:.17g})".format(*p) for p in points)
        root.write_text(text[:match.start(2)]+replacement+text[match.end(2):])
        run(["checkMesh"],case,"log.checkMesh")
    _,_,patches=mesh_lists(case)
    exact=sphere_pressure if sphere else lambda x:plane_pressure(x,omega=2*np.pi*frequency,homogeneous=scenario=="homogeneous")
    if scenario=="sphereWedge": exact=lambda x:sphere_pressure(x[:,[1,0,2]])
    fields={"Pre":("volScalarField","[1 -1 -2 0 0 0 0]","0"),
            "Pim":("volScalarField","[1 -1 -2 0 0 0 0]","0"),
            "pa":("volScalarField","[1 -1 -2 0 0 0 0]","0"),
            "alpha.water":("volScalarField","[0 0 0 0 0 0 0]","0"),
            "alphaf":("surfaceScalarField","[0 0 0 0 0 0 0]","0"),
            "rho":("volScalarField","[1 -3 0 0 0 0 0]","1.2"),
            "sigma":("volTensorField","[0 0 -1 0 0 0 0]","(0 0 0 0 0 0 0 0 0)")}
    for name in ["U","Ure","Uim"]: fields[name]=("volVectorField","[0 1 -1 0 0 0 0]","(0 0 0)")
    for name,(kind,dim,zero) in fields.items():
        body=header(name,kind)+f"dimensions {dim};\ninternalField uniform {zero};\nboundaryField {{\n"
        for patch,(typ,centres) in patches.items():
            if typ in ("empty","symmetryPlane","wedge"):
                entry=f"type {typ}; value uniform {zero};"
            elif name in ("Pre","Pim"):
                pressure=exact(centres)
                entry="type fixedValue; value "+values(pressure.real if name=="Pre" else pressure.imag)+";"
            elif name in ("alphaf","sigma"):
                entry=f"type calculated; value uniform {zero};"
            else: entry=f"type zeroGradient;"
            body+=f"{patch} {{ {entry} }}\n"
        body+="}\n"
        (case/"0"/name).write_text(body)
    if sphere:
        geom=f"type ellipsoid; origin (0 0 0); effectiveRadius .001; horizontalLongAxis {0.001*aspect};"
    else: geom="type plane; origin (0.0001 0 0); normal (1 0 0);"
    (case/"system/setAlphaFieldDict").write_text(header("setAlphaFieldDict")+"field alpha.water;\n"+geom+"\n")
    if scenario!="homogeneous": run(["setAlphaField"],case,"log.setAlphaField")
    run(["postProcess","-func","writeCellCentres","-time","0"],case,"log.centres")
    run(["postProcess","-func","writeCellVolumes","-time","0"],case,"log.volumes")
    return exact

def metrics(case,prepared,exact):
    centres=field(prepared/"0/C")
    V=field(prepared/"0/V",len(centres))
    p=field(case/"1/Pre",len(V))+1j*field(case/"1/Pim",len(V))
    reference=exact(centres); error=abs(p-reference)
    alpha=field(prepared/"0/alpha.water",len(V))
    result={"pressureRelL2":float(np.sqrt(np.sum(V*error**2)/np.sum(V*abs(reference)**2)))}
    for name,weight in [("liquid",alpha),("gas",1-alpha),("interface",((alpha>1e-8)&(alpha<1-1e-8)).astype(float))]:
        denom=np.sum(V*weight*abs(reference)**2)
        result[name+"RelL2"]=float(np.sqrt(np.sum(V*weight*error**2)/denom)) if denom>0 else None
    for name in ["Ure","Uim","pr","momFlux"]: field(case/"1"/name,len(V))
    return p,result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--scenario",choices=["homogeneous","oblique","sphere","sphereWedge","ellipsoid"],required=True)
    ap.add_argument("--sizes",default="32,48,64")
    ap.add_argument("--angles",default="15,45,75")
    ap.add_argument("--ranks",default="1,4")
    ap.add_argument("--mesh-style",choices=["regular","perturbed"],default="perturbed")
    args=ap.parse_args();args.output=args.output.resolve();args.output.mkdir(parents=True,exist_ok=False)
    exe=Path(shutil.which("acousticHelmholtzFoam")).resolve()
    lib=Path(os.environ["FOAM_USER_LIBBIN"])/"libacousticInterface.so"
    hashes=[sha(exe),sha(lib)]
    source_snapshot(args.output,exe,lib)
    results=[]
    for n in map(int,args.sizes.split(",")):
      for angle in (list(map(float,args.angles.split(","))) if args.scenario=="oblique" else [0]):
        prepared=args.output/f"prepared-N{n}-a{angle}"
        exact=prepare(prepared,args.scenario,n,angle,aspect=2 if args.scenario=="ellipsoid" else 1,mesh_style=args.mesh_style)
        input_snapshot(prepared)
        reference={}
        for mode in COMBINATIONS:
          for ranks in map(int,args.ranks.split(",")):
            case=args.output/f"N{n}-a{angle}-{mode}-np{ranks}"
            entry=dict(N=n,angle=angle,mode=mode,ranks=ranks,meshStyle=args.mesh_style)
            try:
                shutil.copytree(prepared,case);configure(case,mode)
                seconds,residuals=solve(case,exe,ranks)
                p,m=metrics(case,prepared,exact)
                if ranks==1: reference[mode]=p
                if mode in reference:
                    m["parallelDifference"]=float(np.linalg.norm(p-reference[mode])/np.linalg.norm(reference[mode]))
                    if mode!="legacy" and m["parallelDifference"]>1e-9: raise RuntimeError("Parallel disagreement")
                if args.scenario=="ellipsoid":
                    m={k:v for k,v in m.items() if not k.endswith("RelL2")}
                entry.update(completed=True,seconds=seconds,residuals=residuals,**m)
            except Exception as e: entry.update(completed=False,error=str(e))
            results.append(entry); print(json.dumps(entry),flush=True)
            (args.output/"results.json").write_text(json.dumps(results,indent=2,allow_nan=False))
    if [sha(exe),sha(lib)]!=hashes: raise RuntimeError("Binary changed during comparison")
    (args.output/"binaries.json").write_text(json.dumps({"solver":hashes[0],"library":hashes[1]},indent=2))
    failures=sum(not x.get("completed",False) for x in results)
    print("ADVANCED_COMPARISONS_COMPLETED",f"failures={failures}",flush=True)
    if failures: raise SystemExit(1)
if __name__=="__main__":main()
