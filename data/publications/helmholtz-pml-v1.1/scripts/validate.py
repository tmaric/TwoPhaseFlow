#!/usr/bin/env python3
"""Validate portable CPC revision data, checksums and manuscript consistency."""
from pathlib import Path
import csv,json,math,hashlib,re,sys
ROOT=Path(__file__).resolve().parents[1]
S=ROOT/'secondary'
checks=[]
def check(ok,message):
 if not ok: raise AssertionError(message)
 checks.append(message)
def rows(path,delimiter=','):
 return list(csv.DictReader((ROOT/path).open(),delimiter=delimiter))
def load(path): return json.loads((ROOT/path).read_text())
def tex(name): return (ROOT/'manuscript/helmholtz'/name).read_text()
def finite(values): return all(math.isfinite(float(x)) for x in values)
def table_rows(name):
 out=[]
 for line in tex(name).splitlines():
  if '&' in line and re.match(r'\s*(?:\d|Orthogonal|Interior-warped|Coarse|Medium|Fine)',line):
   out.append([c.strip().replace('\\\\','').strip() for c in line.split('&')])
 return out
def number(cell):
 cell=cell.strip('$ ').replace('\\,','')
 if '\\times' in cell:
  a,b=cell.split('\\times'); return float(a)*10**int(re.search(r'10\^\{([^}]+)\}',b)[1])
 return float(cell)
def close(a,b,rtol=6e-4,atol=0): return abs(float(a)-float(b))<=max(atol,rtol*abs(float(b)))

def main():
 manifest=load('checksums.json')
 for rel,want in manifest.items():
  p=ROOT/rel
  check(p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==want,'checksum '+rel)
 state=load('hpc-status.json')
 check(state['all_complete'] and len(state['completed'])==10,'all ten campaigns completed')
 check(state['source']=='7f00fc55a67aefbcfa88bdb46a4fc378bac0db0c','immutable software revision')
 binary_hashes=set()
 for task in state['completed']:
  text=(S/task/'binary-sha256.txt').read_text()
  binary_hashes.add(tuple(line.split()[0] for line in text.splitlines()))
 check(len(binary_hashes)==1,'identical solver and library binaries across campaigns')

 homogeneous=rows('secondary/homogeneous/studyResults/convergence/metrics_summary.csv')
 check(len(homogeneous)==24,'24 homogeneous cases')
 for r in table_rows('homogeneous_convergence_table.tex'):
  family='orthogonal' if r[0]=='Orthogonal' else 'warpedInterior'
  h=next(v for v in homogeneous if v['meshFamily']==family and v['boundaryMode']==r[1].lower() and float(v['cellsPerWavelength'])==float(r[2]))
  for i,key in [(3,'pressureRelL2'),(5,'velocityRelL2')]: check(close(h[key],r[i]),f'homogeneous {family} {r[1]} N{r[2]} {key}')
  for i,key in [(4,'pressureOrder'),(6,'velocityOrder')]: check(close(h[key],r[i],atol=.0051),f'homogeneous {family} N{r[2]} {key}')

 layered=load('secondary/layered/comparison/results.json')
 check(len(layered)==66 and all(r['completed'] for r in layered),'66 layered serial/4/8-rank comparisons')
 check(max(r['maxRelativeResidual'] for r in layered)<1e-10,'layered residual gate')
 profiles=rows('secondary/layered/profile-metrics.csv')
 for name in ['pml_sigma_sensitivity_table.tex','pml_mesh_sensitivity_table.tex']:
  for r in table_rows(name):
   is_sigma='sigma' in name
   p=next(v for v in profiles if int(v['N'])==(8000 if is_sigma else int(r[0])) and int(v['sigmaMax'])==(int(r[0]) if is_sigma else 500000) and v['order']=='forward')
   for key,i in [('pressureRelL2',1 if is_sigma else 2),('realRelL2',2 if is_sigma else 4),('imagRelL2',3 if is_sigma else 5)]:
    check(close(p[key],number(r[i])),f'{name} {r[0]} {key}')

 sphere=load('secondary/sphere/studyResults/results.json')
 check(len(sphere)==16 and all(r['completed'] for r in sphere),'16 penetrable-sphere comparisons')
 for r in table_rows('geometric_comparison_table.tex'):
  n=int(r[0])
  for mode,col in [('legacy',1),('geometric',2)]:
   v=next(v for v in sphere if v['N']==n and v['mode']==mode and v['ranks']==1)
   check(close(v['pressureRelL2'],number(r[col]),rtol=7e-5),f'sphere N{n} {mode} pressure')
  v=next(v for v in sphere if v['N']==n and v['mode']=='geometric' and v['ranks']==8)
  check(v['parallelDifference']<2e-14,f'sphere N{n} geometric MPI agreement')
  v=next(v for v in sphere if v['N']==n and v['mode']=='legacy' and v['ranks']==8)
  check(v['parallelDifference']<2e-14 if n==16 else close(v['parallelDifference'],number(r[3]),rtol=1e-4),f'sphere N{n} legacy MPI difference')

 sphere_pml=[]
 for study,count in [('main',24),('soft',16),('ranks',4)]:
  p=ROOT/f'secondary/sphere-pml/studyResults/{study}/results.json'
  if not p.exists():
   alternatives=list((ROOT/'secondary/sphere-pml/studyResults').glob(f'*{study}*/results.json'))
   check(len(alternatives)==1,'sphere-PML study '+study); p=alternatives[0]
  data=json.loads(p.read_text())
  check(len(data)==count and all(r['completed'] for r in data),f'{count} sphere-PML {study} cases')
  sphere_pml.extend(data)
  if study=='main': main_pml=data
 for r in table_rows('sphere_pml_table.tex'):
  for material,c0 in [('waterInAir',1),('oilInWater',3)]:
   for trunc,col in [('dirichlet',c0),('pml',c0+1)]:
    v=next(v for v in main_pml if v['N']==int(r[0]) and v['material']==material and v['truncation']==trunc)
    check(close(100*v['physicalRelL2'],r[col],atol=.00051),f'sphere-PML N{r[0]} {material} {trunc}')
 rank_comparison=rows('secondary/sphere-pml/studyResults/ranks8/rank_comparison.csv')
 check(len(rank_comparison)==4 and max(float(v['relL2Difference']) for v in rank_comparison)<1e-12,'sphere-PML serial/MPI reproduction')
 rigid_difference=rows('secondary/sphere-pml/studyResults/main/rigid_difference.csv')
 for v in rigid_difference:
  check(close(v['relL2Difference'],1.9e-4 if v['material']=='waterInAir' else .14,rtol=.035),'analytical transmitting/rigid reference difference '+v['material'])

 piston=rows('secondary/piston/meshConvergence/metrics_summary.csv')
 check(len(piston)==5,'five piston meshes')
 for r,p in zip(table_rows('piston_convergence_table.tex'),piston):
  check(int(r[0])==int(p['cellsPerWavelength']),'piston resolution')
  for col,key in [(2,'field_relL2'),(4,'axis_relL2'),(5,'axis_relLinf')]: check(close(p[key],number(r[col])),f'piston N{r[0]} {key}')

 gorkov=rows('secondary/gorkov/meshConvergence/mesh_convergence.tsv','\t')
 force_diffs=[]
 check(len(gorkov)==3,'three rigid-sphere mesh levels')
 for r,g in zip(table_rows('gorkov_mesh_convergence_table.tex'),gorkov):
  check(int(r[3])==int(g['n_surface_faces']),'rigid-sphere surface faces '+r[0])
  force=float(g['numerical_force_y_N'])*1e15
  force_diffs.append(abs(force-number(r[4]))/abs(number(r[4])))
  # Cross-build force cancellation changes the last printed digits; document it.
  check(close(force,number(r[4]),rtol=1e-4),'rigid-sphere force within 0.01% '+r[0])
  check(close(100*float(g['difference_from_gorkov']),number(r[5]),atol=.002),'rigid-sphere error within 0.002 percentage points '+r[0])
 radii=rows('secondary/gorkov/gorkovStudy/gorkov_comparison.tsv','\t')
 positions=rows('secondary/gorkov/positionSweep/position_sweep.tsv','\t')
 check(len(radii)==7 and len(positions)==11,'seven radii and eleven positions')
 check(finite(r['numerical_force_y_N'] for r in radii+positions),'finite rigid-sphere forces')
 check(max(float(r['absolute_relative_error']) for r in radii)<.015,'rigid-sphere radius errors below 1.5%')

 height=[]
 for series,key in [('empty','empty'),('1','ab1'),('2','ab2'),('3','ab3')]:
  data=rows(f'secondary/height-{series}/sweepResults/heightSweepResults.csv')
  check(len(data)==100,f'100 height points: {key}')
  for i,r in enumerate(data):
   check(r['series_key']==key and int(r['index'])==i,f'height index {key} {i}')
   check(abs(float(r['height_fac'])-(.46+.18*i/99))<3e-15,f'height grid {key} {i}')
   wedge=float(r['Fn_reflector1_wedge_N']); axis=float(r['Fn_reflector1_axisym_N'])
   check(math.isfinite(axis) and close(axis,360*wedge,rtol=1e-13),f'axisymmetric force factor {key} {i}')
  height.extend(data)
 check(len(height)==400,'400 total height points')
 report={'status':'passed','checks':len(checks),'campaigns':10,'numerical_cases':576,'integration_tests':20,
  'max_rigid_sphere_force_relative_difference_from_printed_workstation_table':max(force_diffs),
  'note':'Printed workstation tables remain primary; HPC reproduction is independent. Force last digits vary across builds. Workstation performance values are transcriptions, not HPC timings.'}
 print(json.dumps(report,indent=2))
 (ROOT/'VALIDATION.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__': main()
