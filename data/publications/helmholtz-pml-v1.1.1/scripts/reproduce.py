#!/usr/bin/env python3
"""Plot the independent HPC reproduction using portable CSV/TSV/JSON only."""
from pathlib import Path
import argparse,csv,json,gzip,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def rows(path,delimiter=','): return list(csv.DictReader((ROOT/path).open(),delimiter=delimiter))
def data(path): return json.loads((ROOT/path).read_text())
def nums(r,k): return np.array([float(v[k]) for v in r])
def save(fig,name,out):
 fig.tight_layout(); fig.savefig(out/(name+'.pdf')); fig.savefig(out/(name+'.png'),dpi=150); plt.close(fig)
def line(ax,r,x,y,label,log=False):
 fn=ax.loglog if log else ax.plot
 fn(nums(r,x),nums(r,y),'o-',label=label,markersize=3)
def main():
 p=argparse.ArgumentParser(); p.add_argument('--output',required=True,type=Path); args=p.parse_args()
 out=args.output.resolve()
 if out.exists(): raise SystemExit('Use a new output directory')
 out.mkdir(parents=True)
 h=rows('secondary/homogeneous/studyResults/convergence/metrics_summary.csv')
 for key in ['pressureRelL2','velocityRelL2']:
  fig,ax=plt.subplots(figsize=(6,4))
  for f in ['orthogonal','warpedInterior']:
   for b in ['dirichlet','mixed']:
    r=[v for v in h if v['meshFamily']==f and v['boundaryMode']==b]
    line(ax,r,'hOverLambda',key,f+' '+b,True)
  ax.set(xlabel='$h/\\lambda$',ylabel=key); ax.legend(fontsize=8); save(fig,'homogeneous-'+key,out)
 profiles=rows('secondary/layered/profile-metrics.csv')
 fig,axes=plt.subplots(2,1,figsize=(6,6))
 for ax,order in zip(axes,['forward','reverse']):
  r=rows(f'secondary/layered/pressure-N8000-s500000-{order}.csv')
  ax.plot(nums(r,'x_m'),nums(r,'exact_absP_Pa'),'k-',label='Reference')
  ax.plot(nums(r,'x_m'),nums(r,'absP_Pa'),'--',label='HPC finite volume')
  ax.set(xlabel='$x$ (m)',ylabel='$|P|$ (Pa)',title=order); ax.legend()
 save(fig,'layered-profiles',out)
 for name,x,r in [('damping','sigmaMax',[v for v in profiles if v['N']=='8000' and v['order']=='forward']),('mesh','h_m',[v for v in profiles if v['N']!='8000'])]:
  fig,ax=plt.subplots(figsize=(6,4))
  for y in ['pressureRelL2','realRelL2','imagRelL2']: line(ax,r,x,y,y,True)
  ax.set(xlabel=x,ylabel='Relative error'); ax.legend(); save(fig,'layered-'+name,out)
 sphere=data('secondary/sphere/studyResults/results.json')
 fig,axes=plt.subplots(1,2,figsize=(10,4))
 for mode in ['legacy','geometric']:
  r=[v for v in sphere if v['mode']==mode and v['ranks']==1]
  line(axes[0],r,'N','pressureRelL2',mode,True)
  r=[v for v in sphere if v['mode']==mode and v['ranks']==8]
  axes[1].loglog(nums(r,'N'),np.maximum(nums(r,'parallelDifference'),1e-16),'o-',label=mode)
 for ax,y in zip(axes,['Pressure error','Serial/eight-rank difference']): ax.set(xlabel='$N$',ylabel=y); ax.legend()
 save(fig,'sphere-face-comparison',out)
 sp=data('secondary/sphere-pml/studyResults/main/results.json')
 fig,axes=plt.subplots(1,2,figsize=(10,4))
 for ax,material in zip(axes,['waterInAir','oilInWater']):
  for trunc in ['dirichlet','pml','pmlCartesian']:
   r=[v for v in sp if v['material']==material and v['truncation']==trunc]
   line(ax,r,'N','physicalRelL2',trunc,True)
  ax.set(xlabel='$N$',ylabel='Pressure error',title=material); ax.legend()
 save(fig,'sphere-pml',out)
 piston=rows('secondary/piston/meshConvergence/metrics_summary.csv')
 fig,ax=plt.subplots(figsize=(6,4))
 for key in ['field_relL2','axis_relL2','axis_relLinf']: line(ax,piston,'h_over_lambda',key,key,True)
 ax.set(xlabel='$h/\\lambda$',ylabel='Relative error'); ax.legend(); save(fig,'piston-convergence',out)
 for name,x,ys in [('onAxisComparison','z_over_rayleigh',['p_sim_over_p0','p_exact_over_p0']),('arcComparison','theta_deg',['SPL_sim_dB','SPL_exact_dB'])]:
  r=rows('secondary/piston/meshConvergence/N80/analyticalCompare/'+name+'.csv')
  fig,ax=plt.subplots(figsize=(6,4))
  for key in ys: ax.plot(nums(r,x),nums(r,key),label=key)
  ax.set(xlabel=x,ylabel='Pressure amplitude' if 'Axis' in name else 'SPL (dB)'); ax.legend(); save(fig,'piston-'+name,out)
 with gzip.open(ROOT/'secondary/piston/meshConvergence/N80/analyticalCompare/fieldComparison.csv.gz','rt') as f: r=list(csv.DictReader(f))
 num=nums(r,'Pre_sim_pa')+1j*nums(r,'Pim_sim_pa'); exact=nums(r,'Pre_exact_pa')+1j*nums(r,'Pim_exact_pa')
 fig,axes=plt.subplots(1,2,figsize=(10,4))
 # A scatter map displays the cell-centre data without introducing interpolation.
 for ax,values,title in zip(axes,[abs(exact),abs(num-exact)],['Exact amplitude (Pa)','Complex-pressure error (Pa)']):
  im=ax.scatter(nums(r,'z_m'),nums(r,'r_m'),c=values,s=1,rasterized=True,cmap='viridis')
  fig.colorbar(im,ax=ax); ax.set(xlabel='$z$ (m)',ylabel='$r$ (m)',title=title)
 save(fig,'piston-field-error',out)
 for name,path,x in [('radii','gorkovStudy/gorkov_comparison.tsv','radius_um'),('positions','positionSweep/position_sweep.tsv','position_mm')]:
  r=rows('secondary/gorkov/'+path,'\t'); fig,ax=plt.subplots(figsize=(6,4))
  for key in ['numerical_force_y_N','gorkov_force_y_N']: ax.plot(nums(r,x),nums(r,key)*1e15,'o-',label=key,markersize=3)
  ax.set(xlabel=x,ylabel='Force (fN)'); ax.legend(fontsize=8); save(fig,'rigid-sphere-'+name,out)
 height=[]
 for s in ['empty','1','2','3']: height+=rows(f'secondary/height-{s}/sweepResults/heightSweepResults.csv')
 fig,ax=plt.subplots(figsize=(6,4))
 for key in ['empty','ab1','ab2','ab3']:
  r=[v for v in height if v['series_key']==key]
  ax.plot(nums(r,'height_mm'),nums(r,'Fn_reflector1_axisym_N'),label=r[0]['series_label'])
 ax.set(xlabel='Gap height (mm)',ylabel='Force on reflector (N)'); ax.legend(fontsize=8); save(fig,'levitator-height-sweep',out)
 # Preserve source tables separately: they contain workstation values, not HPC timings.
 for p in (ROOT/'manuscript/helmholtz').glob('*table.tex'):
  (out/'workstation-reference-tables').mkdir(exist_ok=True); shutil.copy2(p,out/'workstation-reference-tables'/p.name)
 print(f'Generated {len(list(out.glob("*.pdf")))} numerical plots in {out}')
if __name__=='__main__': main()
