from pathlib import Path
import shutil, json, subprocess, sys, csv, hashlib, importlib.util
import numpy as np

root=Path('/work/scratch/tm83tomy/ultrasound-cpc-20261115')
repo=root/'TwoPhaseFlow'
dest=root/'portable-hpc'
dest.mkdir(exist_ok=True)
base=repo/'run/acousticTests/FrequencyDomainTests'
templates={'homogeneous':'homogeneousPlaneWave2D','layered':'layeredInterface1D',
 'sphere':'penetrableSphere','sphere-pml':'penetrableSpherePML','piston':'pistonRadiation',
 'gorkov':'Andrade2019/WedgeLevitatorDropGorkov',
 **{f'height-{x}':'Andrade2019/WedgeLevitatorDropAlpha' for x in ['empty','1','2','3']}}

def copy(src,out):
 out.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(src,out)

completed=[]
for name,template in templates.items():
 task=root/'results'/name
 if not (task/'COMPLETED').exists():
  print('pending',name); continue
 completed.append(name)
 case=task/'case'
 output=dest/'secondary'/name
 output.mkdir(parents=True,exist_ok=True)
 for f in ['binary-sha256.txt','COMPLETED']:
  copy(task/f,output/f)
 for f in ['height-combination-provenance.json','seed.csv','SUPERSEDED']:
  if (task/f).exists(): copy(task/f,output/f)
 # Only actual study products: do not collect copied historical output from nested templates.
 patterns={
  'homogeneous':['studyResults/convergence/*.csv','postProcessing/homogeneousBaseline/1/*'],
  'layered':['comparison/*.json','comparison/*.csv'],
  'sphere':['studyResults/*.json','studyResults/*.csv'],
  'sphere-pml':['studyResults/*/*.json','studyResults/*/*.csv','studyResults/main/*/reference_check.json'],
  'piston':['meshConvergence/*.csv','meshConvergence/N*/analyticalCompare/*'],
  'gorkov':['meshConvergence/*.tsv','gorkovStudy/*.tsv','positionSweep/*.tsv'],
 }.get(name,['sweepResults/heightSweepResults.csv'])
 source=task if name=='layered' else case
 for pattern in patterns:
  for p in source.glob(pattern):
   if p.is_file(): copy(p, output/p.relative_to(source))
 # Keep the exact postprocessing and case parameters as a source snapshot.
 for p in (base/template).rglob('*'):
  if p.is_file() and (p.suffix in ['.py','.sh'] or p.name in ['Allrun','Allclean','caseParams.sh','prepareCase','runMeshConvergence','runGorkovStudy','runPositionSweep']):
   if any(x in ['studyResults','gorkovStudy','positionSweep','meshConvergence','dynamicCode'] for x in p.relative_to(base/template).parts): continue
   copy(p,dest/'scripts/original'/name/p.relative_to(base/template))
 for log in (root/'logs').glob(f'run-{name}-*.log'):
  copy(log,dest/'execution'/log.name)
 if name=='layered':
  sys.path.insert(0,str(case))
  spec=importlib.util.spec_from_file_location('layered_export',case/'run_face_area_comparison.py')
  mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
  layer_metrics=[]
  configs=[(8000,s,'forward') for s in [10000,50000,100000,200000,500000]]+[(8000,500000,'reverse')]+[(n,500000,'forward') for n in [560,1120,2240,4480,8960]]
  for n,sigma,order in configs:
   tag=f'N{n}-s{sigma}-{order}'; prepared=task/'comparison'/f'prepared-{tag}'
   props=mod.parse_keyvals(prepared/'studyProperties')
   x=props['X_MIN']+(np.arange(n)+.5)*(props['X_MAX']-props['X_MIN'])/n
   num=task/'comparison'/f'{tag}-geometric-np1'/'1'
   p=mod.field(num/'Pre',n)+1j*mod.field(num/'Pim',n); ref=mod.exact(prepared,n)
   np.savetxt(output/f'pressure-{tag}.csv', np.column_stack([x,p.real,p.imag,abs(p),ref.real,ref.imag,abs(ref)]),delimiter=',',header='x_m,Pre_Pa,Pim_Pa,absP_Pa,exact_Pre_Pa,exact_Pim_Pa,exact_absP_Pa',comments='')
   layer_metrics.append(dict(N=n,sigmaMax=sigma,order=order,h_m=(props['X_MAX']-props['X_MIN'])/n,
    pressureRelL2=float(np.linalg.norm(p-ref)/np.linalg.norm(ref)),
    realRelL2=float(np.linalg.norm(p.real-ref.real)/np.linalg.norm(ref.real)),
    imagRelL2=float(np.linalg.norm(p.imag-ref.imag)/np.linalg.norm(ref.imag))))
  with (output/'profile-metrics.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,list(layer_metrics[0])); w.writeheader(); w.writerows(layer_metrics)
 if name=='gorkov':
  for study,pattern in [('meshConvergence','*'),('gorkovStudy','r*um'),('positionSweep','p*mm')]:
   for run in (case/study).glob(pattern):
    f=run/'postProcessing/dropRadiationForce/force.dat'
    if f.is_file(): copy(f,output/study/run.name/'force.dat')

for p in (root/'scripts').iterdir():
 if p.is_file(): copy(p,dest/'execution'/p.name)
for name in ['campaign-jobs.tsv','height-shard-plan.json','height-shard-jobs.tsv','heightSweepResults-combined.csv','environment.txt','build-environment.txt','job-ids.txt']:
 if (root/name).exists(): copy(root/name,dest/'execution'/name)
for p in (root/'logs').glob('build-*.log'): copy(p,dest/'execution'/p.name)
for p in (root/'logs').glob('height-shard-*.log'): copy(p,dest/'execution'/p.name)
state={'source':subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD']).decode().strip(),
 'date':'2026-11-15','completed':completed,'expected':list(templates),'all_complete':len(completed)==10}
(dest/'hpc-status.json').write_text(json.dumps(state,indent=2)+'\n')
print(json.dumps(state,indent=2))
