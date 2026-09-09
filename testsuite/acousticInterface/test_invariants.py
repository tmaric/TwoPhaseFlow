import os,re,shutil,sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parent))
from compare import configure,field,prepare_layered,run,solve

@pytest.fixture(scope="module")
def prepared(tmp_path_factory):
    if not shutil.which("acousticHelmholtzFoam"): pytest.skip("Source the OpenFOAM test environment")
    path=tmp_path_factory.mktemp("prepared")/"case"
    prepare_layered(path,32,.25)
    return path

def clone(prepared,tmp_path):
    case=tmp_path/"case";shutil.copytree(prepared,case);return case

def test_compiled_local_geometry_and_transmission(prepared):
    run(["testAcousticInterface"],prepared,"log.unit")
    assert "ACOUSTIC_INTERFACE_TESTS_PASSED" in (prepared/"log.unit").read_text()

@pytest.mark.parametrize("settings,expected",[
    ("{ areaFraction legacy; flux plicTransmission; }","requires areaFraction plicAverage"),
    ("{ areaFraction typo; }","Unknown acousticInterface"),
    ("{ flux typo; }","Unknown acousticInterface"),
    ("{ plicTransmissionCoeffs { maxStencilRings 0; } }","Invalid plicTransmissionCoeffs"),
    ("{ plicTransmissionCoeffs { svdRelativeTolerance 0; } }","Invalid plicTransmissionCoeffs")])
def test_invalid_configuration(prepared,tmp_path,settings,expected):
    case=clone(prepared,tmp_path)
    run(["foamDictionary","system/fvSchemes","-entry","acousticInterface","-set",settings],case,"log.configure")
    run(["acousticHelmholtzFoam"],case,"log.solver",expect_failure=True)
    assert expected in (case/"log.solver").read_text()

def test_missing_dictionary_is_legacy(prepared,tmp_path):
    case=clone(prepared,tmp_path);solve(case,"acousticHelmholtzFoam",1)
    original=field(case/"1/Pre",32)+1j*field(case/"1/Pim",32)
    configure(case,"legacy");solve(case,"acousticHelmholtzFoam",1)
    explicit=field(case/"1/Pre",32)+1j*field(case/"1/Pim",32)
    assert np.linalg.norm(original-explicit)/np.linalg.norm(original)<1e-12

@pytest.mark.parametrize("bad",["nan","inf","-inf"])
def test_nonfinite_postprocessing_rejected(tmp_path,bad):
    p=tmp_path/"Pre"
    p.write_text(f"internalField nonuniform List<scalar> 3 (1 {bad} 2);")
    with pytest.raises(ValueError,match="Nonfinite"):field(p,3)

def test_wrong_field_length_rejected(tmp_path):
    p=tmp_path/"Pre";p.write_text("internalField nonuniform List<scalar> 3 (1 2);")
    with pytest.raises(ValueError):field(p,3)

def test_phi_sign_independence(prepared,tmp_path):
    values=[]
    for sign in [1,-1]:
        case=tmp_path/str(sign);shutil.copytree(prepared,case);configure(case,"geometry")
        (case/"0/phi").write_text(f"""FoamFile {{ version 2.0; format ascii; class surfaceScalarField; object phi; }}
dimensions [0 3 -1 0 0 0 0]; internalField uniform {sign};
boundaryField {{
 left {{ type calculated; value uniform {sign}; }}
 right {{ type calculated; value uniform {sign}; }}
 yMin {{ type empty; }} yMax {{ type empty; }}
 zMin {{ type empty; }} zMax {{ type empty; }}
}}
""")
        solve(case,"acousticHelmholtzFoam",1);values.append(field(case/"1/alphaf",31))
    np.testing.assert_allclose(values[0],values[1],atol=1e-12,rtol=0)

def test_interface_pml_overlap_rejected(prepared,tmp_path):
    case=clone(prepared,tmp_path);configure(case,"transmission")
    text=(case/"0/sigma").read_text()
    text=re.sub(r"internalField\s+uniform\s+\([^;]+\);",
                "internalField uniform (1 0 0 0 0 0 0 0 0);",text)
    (case/"0/sigma").write_text(text)
    run(["acousticHelmholtzFoam"],case,"log.solver",expect_failure=True)
    assert "Interface/PML overlap" in (case/"log.solver").read_text()

def read_petsc_matrix(path):
    from scipy.sparse import csr_matrix
    with Path(path).open("rb") as f:
        header=np.fromfile(f,dtype=">i4",count=4)
        assert header[0]==1211216,header
        _,rows,cols,nz=header
        lengths=np.fromfile(f,dtype=">i4",count=rows)
        indices=np.fromfile(f,dtype=">i4",count=nz).astype(np.int64)
        data=np.fromfile(f,dtype=">f8",count=nz)
        assert len(data)==nz and np.isfinite(data).all()
    indptr=np.r_[0,np.cumsum(lengths)]
    return csr_matrix((data,indices,indptr),shape=(rows,cols))

def test_actual_petsc_1d_matrix(prepared,tmp_path,monkeypatch):
    case=clone(prepared,tmp_path);configure(case,"transmission")
    monkeypatch.setenv("PETSC_OPTIONS","-ksp_view_mat binary:matrix.dat")
    solve(case,"acousticHelmholtzFoam",1)
    matrix=read_petsc_matrix(case/"matrix.dat").toarray()[1::2,1::2]
    n=32;h=.14/n;area=4e-6
    alpha=field(case/"0/alpha.water",n);rho=alpha*1000+(1-alpha)*1.2
    kappa=alpha/(1000*1500**2)+(1-alpha)/(1.2*343**2)
    diffusion=-matrix/rho[:,None]+np.diag(area*h*(2*np.pi*20000)**2*kappa)
    # Build the same operator independently from geometric path resistances.
    x=(np.arange(n)+.5)*h;xi=.07+.25*h
    liquid=np.maximum(0,x[1:]-np.maximum(x[:-1],xi))
    gas=h-liquid
    g=area/(gas*1.2+liquid*1000)
    expected=np.diag(np.r_[2*area/(1.2*h),g]+np.r_[g,2*area/(1000*h)])
    expected-=np.diag(g,1)+np.diag(g,-1)
    assert np.linalg.norm(diffusion-expected)/np.linalg.norm(expected)<1e-10
    assert np.linalg.eigvalsh(diffusion).min()>0
    neumann=diffusion.copy()
    neumann[0,0]-=2*area/(1.2*h);neumann[-1,-1]-=2*area/(1000*h)
    scale=np.linalg.norm(neumann,2)
    eig=np.linalg.eigvalsh(neumann)
    assert abs(eig[0])<1e-10*scale and eig[1]>1e-10*scale
    assert np.linalg.norm(neumann@np.ones(n))<1e-10*scale
    assert np.linalg.norm(neumann@((-1.)**np.arange(n)))>1e-3*scale
    rng=np.random.default_rng(46);p=rng.normal(size=n)
    assert np.linalg.norm(diffusion@p-expected@p)<1e-10*scale*np.linalg.norm(p)


def test_analytical_sphere_transmission_conditions():
    from advanced import sphere_pressure
    a=.001;h=a*1e-5
    angles=np.linspace(.1,3.0,13)
    normals=np.column_stack([np.cos(angles),np.sin(angles),np.zeros_like(angles)])
    p0=sphere_pressure(a*normals)
    pin=sphere_pressure((a-h)*normals)
    pin2=sphere_pressure((a-2*h)*normals)
    pout=sphere_pressure((a+h)*normals)
    pout2=sphere_pressure((a+2*h)*normals)
    qin=(3*p0-4*pin+pin2)/(2*h*1000)
    qout=(-3*p0+4*pout-pout2)/(2*h*1.2)
    assert np.linalg.norm(qin-qout)/np.linalg.norm(qout)<1e-5
    assert np.linalg.norm(pin-pout)/np.linalg.norm(p0)<1e-5


def test_empty_samples_rejected(tmp_path):
    p=tmp_path/"Pre";p.write_text("internalField nonuniform List<scalar> 0 (); ")
    with pytest.raises(ValueError,match="No field samples"):field(p)


@pytest.mark.parametrize("mode",["geometry","transmission"])
def test_cyclic_patch_rejected(prepared,tmp_path,mode):
    case=clone(prepared,tmp_path)
    p=case/"system/blockMeshDict";text=p.read_text()
    for patch,other in [("left","right"),("right","left")]:
        text=re.sub(r"("+patch+r"\s*\{\s*)type patch;",
                    r"\1type cyclic; neighbourPatch "+other+";",text)
    p.write_text(text)
    run(["blockMesh"],case,"log.cyclicMesh")
    for p in (case/"0").iterdir():
        if not p.is_file(): continue
        text=p.read_text()
        for patch in ["left","right"]:
            text=re.sub(r"("+patch+r"\s*\{).*?\}",r"\1 type cyclic; }",text,flags=re.S)
        p.write_text(text)
    configure(case,mode)
    run(["acousticHelmholtzFoam"],case,"log.solver",expect_failure=True)
    assert "New acoustic interface methods do not support patch" in (case/"log.solver").read_text()


@pytest.mark.parametrize("phase",[0,1])
def test_homogeneous_phase_with_nonzero_pml(prepared,tmp_path,phase):
    results=[]
    for mode in ["legacy","geometry","transmission"]:
        case=tmp_path/mode;shutil.copytree(prepared,case);configure(case,mode)
        run(["foamDictionary","0/alpha.water","-entry","internalField","-set",f"uniform {phase}"],case,"log.pure")
        run(["foamDictionary","0/sigma","-entry","internalField","-set", "uniform (10000 0 0 0 0 0 0 0 0)"],case,"log.pml")
        solve(case,"acousticHelmholtzFoam",1)
        results.append(field(case/"1/Pre",32)+1j*field(case/"1/Pim",32))
    for result in results[1:]:
        assert np.linalg.norm(result-results[0])/np.linalg.norm(results[0])<1e-12
