#!/usr/bin/env python3
"""Validate the portable manuscript secondary-data deposit, without OpenFOAM."""
from pathlib import Path
import argparse, csv, hashlib, json, math, re
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def rows(path):
    with (ROOT/path).open(newline='') as f: return list(csv.DictReader(f))

def col(rr, name): return np.array([float(r[name]) for r in rr])

def near(a,b,rtol=1e-6,atol=0):
    if not np.allclose(a,b,rtol=rtol,atol=atol):
        raise AssertionError(f'Numeric mismatch: max absolute difference {np.max(np.abs(np.asarray(a)-np.asarray(b)))}')

def inside(path):
    p=ROOT/path
    if Path(path).is_absolute() or not p.resolve().is_relative_to(ROOT):
        raise AssertionError(f'Nonportable path: {path}')
    if not p.is_file(): raise AssertionError(f'Missing file: {path}')
    return p

def run(checksums=True):
    checks=[]
    if checksums:
        m=json.loads((ROOT/'checksums.json').read_text())
        actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='checksums.json'}
        assert actual==set(m['files']),('Unexpected/missing payload files',actual.symmetric_difference(m['files']))
        for name,record in m['files'].items():
            p=inside(name)
            assert p.stat().st_size==record['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==record['sha256'],name
        checks.append(f"SHA-256 and sizes checked for {len(actual)} files")
    d=json.loads((ROOT/'data-dictionary.json').read_text())
    files=list((ROOT/'secondary').rglob('*.csv'))
    assert {p.relative_to(ROOT).as_posix() for p in files}==set(d['files'])
    for p in files:
        key=p.relative_to(ROOT).as_posix();record=d['files'][key];rr=rows(key)
        assert len(rr)==record['rows'] and list(rr[0])==record['columns'],key
        for row in rr:
            for name,value in row.items():
                spec=d['columns'][name]
                if value=='':
                    assert name in record['allowed_empty_columns'],(key,name)
                elif spec['type']=='number':
                    assert math.isfinite(float(value)),(key,name,value)
                elif spec['type']=='path': inside(value)
    checks.append(f"All {len(files)} CSV files match the schema, row counts, missing-value rules and finite-value checks")
    prov=json.loads((ROOT/'provenance.json').read_text());items=prov['items']
    assert sorted(x['number'] for x in items if x['kind']=='Figure')==list(range(1,20))
    assert sorted(x['number'] for x in items if x['kind']=='Table')==list(range(1,8))
    assert len({x['label'] for x in items})==26
    manuscript=inside(prov['manuscript_source'])
    assert hashlib.sha256(manuscript.read_bytes()).hexdigest()==prov['manuscript_source_sha256']
    all_text=manuscript.read_text()
    for name in re.findall(r'\\paperinput\{([^}]+)\}',all_text):all_text+='\n'+inside('manuscript/'+name).read_text()
    for item in items:
        assert '\\label{'+item['label']+'}' in all_text,item['label']
        for path in item['data_paths']+item['asset_paths']+[item['reproduce']]:inside(path)
    for path in prov['text_only_support']:inside(path)
    checks.append('All 19 figures, 7 tables and text-only supporting results resolve within this deposit')
    # The figures and error tables share these exact exported sample arrays.
    for r in rows('secondary/layered/metrics.csv'):
        rr=rows(r['profile_path']);assert len(rr)==1200
        x=col(rr,'x_m');assert np.all(np.diff(x)>0);near([x[0],x[-1]],[0,.35],atol=2e-8)
        p=col(rr,'Pre_sim')+1j*col(rr,'Pim_sim');q=col(rr,'Pre_analytic')+1j*col(rr,'Pim_analytic')
        near(np.linalg.norm(p-q)/np.linalg.norm(q),float(r['P_relL2']),rtol=3e-6)
        for name,component in [('Pre',np.real),('Pim',np.imag)]:
            near(np.linalg.norm(component(p-q))/np.linalg.norm(component(q)),float(r[name+'_relL2']),rtol=3e-6)
        near(np.abs(p),col(rr,'p_abs_sim'),rtol=2e-6,atol=1e-8)
    checks.append('All 11 layered pressure norms independently recomputed from the 1200-point profiles')
    comp=rows('secondary/layered/method-comparisons.csv');assert len(comp)==76
    assert np.max(col(comp,'baselineDifference'))<1.95e-12
    for r in comp:
        if r['mode']=='submission':assert r['max_relative_residual'] and float(r['max_relative_residual'])<=1e-10
        else:assert r['mode']=='baseline' and r['max_relative_residual']==''
    for r in rows('secondary/piston/metrics.csv'):
        n=int(r['cellsPerWavelength']);on=rows(f'secondary/piston/on-axis-N{n}.csv');ff=rows(f'secondary/piston/far-field-N{n}.csv')
        assert len(on)==400 and len(ff)==361
        sim=col(on,'p_sim_over_p0');exact=col(on,'p_analytic_over_p0');err=sim-exact
        near(np.linalg.norm(err)/np.linalg.norm(exact),float(r['relL2']),rtol=4e-6)
        near(np.max(np.abs(err))/np.max(np.abs(exact)),float(r['relLinf']),rtol=4e-6)
        sim=col(ff,'p_far_abs_sim_pa');exact=col(ff,'p_far_abs_analytic_pa')
        near(np.linalg.norm(sim-exact)/np.linalg.norm(exact),float(r['farField_pressureMagnitude_relL2']),rtol=1e-6)
        near(col(ff,'SPL_sim_dB'),10*np.log10(.5*(sim/20e-6)**2),rtol=1e-12)
        near(col(ff,'directivity_sim'),sim/sim[0],rtol=1e-12)
    checks.append('Both piston error tables, SPL and directivity recomputed from all five exported profile pairs')
    for file in ['convergence','boundary-skew-diagnostic']:
        rr=rows(f'secondary/homogeneous/{file}.csv')
        for family in set(x['meshFamily'] for x in rr):
            for bc in set(x['boundaryMode'] for x in rr):
                group=sorted([x for x in rr if x['meshFamily']==family and x['boundaryMode']==bc],key=lambda x:float(x['cellsPerWavelength']))
                assert group[0]['pressureOrder']==group[0]['velocityOrder']==''
                for before,after in zip(group,group[1:]):
                    for error,order in [('pressureRelL2','pressureOrder'),('velocityRelL2','velocityOrder')]:
                        q=math.log(float(before[error])/float(after[error]))/math.log(float(before['hOverLambda'])/float(after['hOverLambda']))
                        near(q,float(after[order]),rtol=1e-6)
    checks.append('Homogeneous pressure and velocity orders recomputed for both retained mesh studies')
    k=2*math.pi*25250/343;energy=1/(4*1.2*343**2)
    def force(a,y):return 4*math.pi*np.asarray(a)**3*k*energy*(5/6)*np.sin(2*k*np.asarray(y))
    for file in ['mesh-convergence','radius-sweep','position-sweep']:
        rr=rows(f'secondary/gorkov/{file}.csv');y=col(rr,'position_m') if file=='position-sweep' else .0040003;a=col(rr,'radius_m') if file=='radius-sweep' else 60e-6
        exact=force(a,y);near(exact,col(rr,'gorkov_force_y_N'),rtol=1e-8,atol=1e-25)
        num=col(rr,'numerical_force_y_N')
        if file=='position-sweep':near(np.abs(num-exact)/abs(4*math.pi*a**3*k*energy*(5/6)),col(rr,'difference_normalized_by_peak_force'),rtol=2e-6)
        else:near(np.abs(num-exact)/np.abs(exact),col(rr,'absolute_relative_error' if file=='radius-sweep' else 'difference_from_gorkov'),rtol=2e-6)
    curve=rows('secondary/gorkov/position-analytical-curve.csv');near(force(60e-6,col(curve,'position_m')),col(curve,'gorkov_force_y_N'),atol=1e-29)
    checks.append('Gorkov analytical forces, error measures and the 500-point analytical curve independently checked')
    # Check every displayed table entry against the secondary data at printed precision.
    def table_rows(name):
        t=inside(name).read_text()
        return [l for l in t.splitlines() if '&' in l and l.strip().endswith(r'\\') and re.match(r'\s*(\d|Orthogonal|Interior-warped|Coarse|Medium|Fine)',l)]
    def fields(line):return [v.strip().replace(r'\\','').strip() for v in line.split('&')]
    for item in [x for x in items if x['kind']=='Table']:
        n=item['number'];lines=table_rows(item['asset_paths'][0])
        if n==1:
            assert len(lines)==16;rr=rows(item['data_paths'][0])
            for line in lines:
                c=fields(line);r=next(r for r in rr if r['meshFamily']==('orthogonal' if c[0]=='Orthogonal' else 'warpedInterior') and r['boundaryMode']==c[1].lower() and int(float(r['cellsPerWavelength']))==int(c[2]))
                for i,k in [(3,'pressureRelL2'),(5,'velocityRelL2')]:near(float(c[i]),float(r[k]),rtol=5e-4)
                for i,k in [(4,'pressureOrder'),(6,'velocityOrder')]:near(float(c[i]),float(r[k]),rtol=0,atol=.005)
        elif n in [2,3]:
            rr=rows('secondary/layered/metrics.csv');rr=sorted([r for r in rr if (int(r['N'])==8000 and r['ordering']=='forward') if n==2],key=lambda r:int(r['sigma_max_s_inv'])) if n==2 else sorted([r for r in rr if int(r['N'])!=8000],key=lambda r:int(r['N']))
            assert len(lines)==len(rr)==5
            for line,r in zip(lines,rr):
                c=fields(line)
                assert int(c[0])==int(r['sigma_max_s_inv'] if n==2 else r['N'])
                if n==3:
                    mantissa,exponent=re.fullmatch(r'\$([\d.]+)\\times\s*10\^\{(-?\d+)\}\$',c[1]).groups()
                    near(float(mantissa)*10**int(exponent),.35/int(r['N']),rtol=5e-4)
                vals=[float(v)*10**int(e) for v,e in re.findall(r'([\d.]+)\\times\s*10\^\{(-?\d+)\}',line)][-3:]
                near(vals,[float(r[k]) for k in ['P_relL2','Pre_relL2','Pim_relL2']],rtol=5e-4)
        elif n==4:
            rr=rows('secondary/sphere/area-comparison.csv');assert len(rr)==16 and len(lines)==4
            for line,N in zip(lines,[16,24,32,48]):
                assert int(fields(line)[0])==N
                get=lambda m,r:next(x for x in rr if int(x['N'])==N and x['mode']==m and int(x['ranks'])==r)
                expected=[float(get(m,1)['pressureRelL2']) for m in ['legacy','geometry']]+[float(get(m,8)['parallelDifference']) for m in ['legacy','geometry']]
                vals=[float(v)*10**int(e) for v,e in re.findall(r'([\d.]+)\\times\s*10\^\{(-?\d+)\}',line)];near(vals,expected,rtol=5e-5)
        elif n in [5,6]:
            rr=rows('secondary/piston/metrics.csv');assert len(lines)==5
            for line,r in zip(lines,rr):
                c=fields(line);assert int(c[0])==int(r['cellsPerWavelength'])
                near(float(c[1]),float(r['h_over_lambda']),rtol=0,atol=.00005)
                for i,k in enumerate(['relL2','relLinf'] if n==5 else ['farField_pressureMagnitude_relL2'],2):near(float(c[i]),float(r[k]),rtol=5e-4)
        elif n==7:
            rr=rows('secondary/gorkov/mesh-convergence.csv');assert len(lines)==3
            for line,r in zip(lines,rr):
                c=fields(line);assert c[0].lower()==r['level'].lower()
                near(float(c[1]),float(r['h_over_a']),rtol=0,atol=.000005)
                assert int(c[2])==int(r['segments']) and int(c[3])==int(r['n_surface_faces'])
                near(float(c[4]),1e15*float(r['numerical_force_y_N']),rtol=0,atol=.0000006)
                near(float(c[5]),100*float(r['difference_from_gorkov']),rtol=0,atol=.0006)
    checks.append('Every numerical row of all seven manuscript tables matches its portable secondary-data source')
    return dict(status='passed',checks=checks,figures=19,tables=7,csv_files=len(files))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--skip-checksums',action='store_true');args=p.parse_args()
    print(json.dumps(run(not args.skip_checksums),indent=2))
