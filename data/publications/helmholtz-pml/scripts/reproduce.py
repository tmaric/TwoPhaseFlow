#!/usr/bin/env python3
"""Replot every numerical manuscript figure and regenerate every table from CSV."""
from pathlib import Path
import argparse, csv, json, sys
sys.dont_write_bytecode=True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from validate import ROOT, rows, col, run

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    out=args.output.resolve()
    if out.is_relative_to(ROOT):raise ValueError('Use an output directory outside the immutable deposit.')
    report=run();out.mkdir(parents=True,exist_ok=False)
    def save(fig,n):
        fig.tight_layout();fig.savefig(out/f'figure{n:02d}.png',dpi=200);plt.close(fig)
    rr=rows('secondary/homogeneous/convergence.csv')
    fig,axes=plt.subplots(1,2,figsize=(9,3.8))
    for family in ['orthogonal','warpedInterior']:
        for bc in ['dirichlet','mixed']:
            group=[r for r in rr if r['meshFamily']==family and r['boundaryMode']==bc]
            for ax,key in zip(axes,['pressureRelL2','velocityRelL2']):
                ax.loglog(col(group,'hOverLambda'),col(group,key),'o-' if bc=='dirichlet' else 's--',label=f'{family}, {bc}')
                ax.set_xlabel(r'$h/\lambda$');ax.set_ylabel(key);ax.grid(True,which='both',alpha=.25)
    for ax in axes:ax.invert_xaxis()
    axes[0].legend(fontsize=7);save(fig,5)
    forward=rows('secondary/layered/profiles/N8000-forward-sigma500000.csv');reverse=rows('secondary/layered/profiles/N8000-reverse-sigma500000.csv')
    def profile(ax,rr,kind):
        x=col(rr,'x_m');sim='p_abs_sim' if kind=='abs' else kind+'_sim';ana='p_abs_analytic' if kind=='abs' else kind+'_analytic'
        ax.plot(x,col(rr,ana),'k-',label='Analytical');ax.plot(x,col(rr,sim),'--',label='Numerical')
        ax.axvspan(.2,.35,color='0.9');ax.set_xlim(0,.35);ax.set_ylabel((r'$|P|$' if kind=='abs' else kind)+' [Pa]');ax.grid(alpha=.25);ax.legend(fontsize=8)
    fig,axes=plt.subplots(2,1,sharex=True,figsize=(7,5.4))
    for ax,kind in zip(axes,['Pre','Pim']):profile(ax,forward,kind)
    axes[-1].set_xlabel('x [m]');save(fig,7)
    fig,ax=plt.subplots(figsize=(7,3.7));profile(ax,forward,'abs');ax.set_xlabel('x [m]');save(fig,8)
    fig,axes=plt.subplots(3,1,sharex=True,figsize=(7,7.6))
    for ax,kind in zip(axes,['Pre','Pim','abs']):profile(ax,reverse,kind)
    axes[-1].set_xlabel('x [m]');save(fig,9)
    lm=rows('secondary/layered/metrics.csv')
    for n,group,xkey in [(10,sorted([r for r in lm if int(r['N'])==8000 and r['ordering']=='forward'],key=lambda r:int(r['sigma_max_s_inv'])),'sigma_max_s_inv'),(11,sorted([r for r in lm if int(r['N'])!=8000],key=lambda r:int(r['N'])),'N')]:
        fig,ax=plt.subplots(figsize=(7,4.3))
        for key,mark in [('P_relL2','o'),('Pre_relL2','s'),('Pim_relL2','^')]:ax.loglog(col(group,xkey),col(group,key),marker=mark,label=key)
        ax.set_xlabel(r'$\sigma_{\max}$ [s$^{-1}$]' if n==10 else r'$N_x$');ax.set_ylabel('Relative pressure error');ax.grid(True,which='both',alpha=.25);ax.legend();save(fig,n)
    fig,ax=plt.subplots(figsize=(7,4.3))
    for n in [20,30,40,60,80]:
        rr=rows(f'secondary/piston/on-axis-N{n}.csv')
        if n==20:ax.plot(col(rr,'z_over_rayleigh'),col(rr,'p_analytic_over_p0'),'k-',label='Analytical')
        ax.plot(col(rr,'z_over_rayleigh'),col(rr,'p_sim_over_p0'),'--',label=f'N={n}')
    ax.set_xscale('log');ax.set_xlim(5e-4,1);ax.set_xlabel(r'$z/R_0$');ax.set_ylabel(r'$|P|/(\rho c u_0)$');ax.grid(alpha=.25);ax.legend();save(fig,13)
    fig,ax=plt.subplots(figsize=(7,4.3))
    for n in [20,30,40,60,80]:
        rr=rows(f'secondary/piston/far-field-N{n}.csv')
        if n==20:ax.plot(col(rr,'theta_deg'),col(rr,'SPL_analytic_dB'),'k-',label='Analytical')
        ax.plot(col(rr,'theta_deg'),col(rr,'SPL_sim_dB'),'--',label=f'N={n}')
    ax.set_xlim(0,90);ax.set_xlabel('Angle [degrees]');ax.set_ylabel(r'SPL [dB re 20 $\mu$Pa]');ax.grid(alpha=.25);ax.legend();save(fig,14)
    rr=rows('secondary/piston/far-field-N80.csv');theta=np.deg2rad(col(rr,'theta_deg'));theta=np.concatenate((-theta[:0:-1],theta))
    fig,ax=plt.subplots(figsize=(6,5.2),subplot_kw={'projection':'polar'});ax.set_theta_zero_location('N');ax.set_theta_direction(-1)
    for key,label,style in [('SPL_analytic_dB','Analytical','-'),('SPL_sim_dB','Numerical','--')]:
        y=col(rr,key);y=np.concatenate((y[:0:-1],y));ax.plot(theta,y-y.max(),style,label=label)
    ax.set_thetamin(-90);ax.set_thetamax(90);ax.set_rlim(-40,0);ax.legend();save(fig,15)
    rr=rows('secondary/gorkov/mesh-convergence.csv');fig,axes=plt.subplots(1,2,figsize=(7.2,3.2));x=col(rr,'h_over_a')
    axes[0].axhline(float(rr[0]['gorkov_force_y_N'])*1e15,color='black',label='Gorkov');axes[0].plot(x,col(rr,'numerical_force_y_N')*1e15,'o-',label='Numerical');axes[0].set_ylabel('Force [fN]');axes[0].legend()
    axes[1].semilogy(x,col(rr,'difference_from_gorkov')*100,'o-');axes[1].set_ylabel('Relative difference [%]')
    for ax in axes:ax.invert_xaxis();ax.set_xlabel(r'$h_s/a$');ax.grid(alpha=.25)
    save(fig,17)
    rr=rows('secondary/gorkov/radius-sweep.csv');fig,axes=plt.subplots(2,1,sharex=True,figsize=(6.4,6.2),gridspec_kw={'height_ratios':(2,1)})
    axes[0].plot(col(rr,'ka'),col(rr,'gorkov_force_y_N')*1e15,'k-',label='Gorkov');axes[0].plot(col(rr,'ka'),col(rr,'numerical_force_y_N')*1e15,'o',label='Numerical');axes[0].set_ylabel('Force [fN]');axes[0].legend()
    axes[1].plot(col(rr,'ka'),col(rr,'absolute_relative_error')*100,'o-');axes[1].set_ylabel('Relative error [%]');axes[1].set_xlabel(r'$ka$')
    for ax in axes:ax.grid(alpha=.25)
    save(fig,18)
    rr=rows('secondary/gorkov/position-sweep.csv');curve=rows('secondary/gorkov/position-analytical-curve.csv');fig,axes=plt.subplots(2,1,sharex=True,figsize=(6.4,6.2),gridspec_kw={'height_ratios':(2,1)})
    axes[0].plot(col(curve,'position_m')*1000,col(curve,'gorkov_force_y_N')*1e15,'k-',label='Gorkov');axes[0].plot(col(rr,'position_m')*1000,col(rr,'numerical_force_y_N')*1e15,'o',label='Numerical');axes[0].set_ylabel('Force [fN]');axes[0].legend()
    axes[1].plot(col(rr,'position_m')*1000,col(rr,'difference_normalized_by_peak_force')*100,'o-');axes[1].set_ylabel('Peak-normalized error [%]');axes[1].set_xlabel('Sphere centre y [mm]')
    for ax in axes:ax.grid(alpha=.25)
    save(fig,19)
    # Regenerate LaTeX tables numerically, retaining useful column names and units.
    def latex(n,headers,values):
        text=['\\begin{tabular}{'+'l'*len(headers)+'}', '\\hline', ' & '.join(headers)+r' \\', '\\hline']
        text+=[' & '.join(str(v) for v in row)+r' \\' for row in values];text+=['\\hline','\\end{tabular}']
        (out/f'table{n:02d}.tex').write_text('\n'.join(text)+'\n')
    rr=[r for r in rows('secondary/homogeneous/convergence.csv') if float(r['cellsPerWavelength'])>=32]
    latex(1,['Mesh','BC',r'$N_\lambda$',r'$E_P$',r'$q_P$',r'$E_u$',r'$q_u$'],[[r['meshFamily'],r['boundaryMode'],int(float(r['cellsPerWavelength'])),f"{float(r['pressureRelL2']):.3e}",f"{float(r['pressureOrder']):.2f}",f"{float(r['velocityRelL2']):.3e}",f"{float(r['velocityOrder']):.2f}"] for r in rr])
    rr=sorted([r for r in lm if int(r['N'])==8000 and r['ordering']=='forward'],key=lambda r:int(r['sigma_max_s_inv']))
    latex(2,[r'$\sigma_{\max}$ [s$^{-1}$]',r'$E_P$',r'$E_{Pre}$',r'$E_{Pim}$'],[[r['sigma_max_s_inv']]+[f'{float(r[k]):.3e}' for k in ['P_relL2','Pre_relL2','Pim_relL2']] for r in rr])
    rr=sorted([r for r in lm if int(r['N'])!=8000],key=lambda r:int(r['N']))
    latex(3,[r'$N_x$',r'$\Delta x$ [m]',r'$E_P$',r'$E_{Pre}$',r'$E_{Pim}$'],[[r['N'],f"{.35/int(r['N']):.3e}"]+[f'{float(r[k]):.3e}' for k in ['P_relL2','Pre_relL2','Pim_relL2']] for r in rr])
    rr=rows('secondary/sphere/area-comparison.csv');vv=[]
    for N in [16,24,32,48]:
        get=lambda m,p:next(r for r in rr if int(r['N'])==N and r['mode']==m and int(r['ranks'])==p)
        vv.append([N]+[f"{float(get(m,1)['pressureRelL2']):.4e}" for m in ['legacy','geometry']]+[f"{float(get(m,8)['parallelDifference']):.4e}" for m in ['legacy','geometry']])
    latex(4,['N','Original $E_P$','Geometric $E_P$','Original $D_8$','Geometric $D_8$'],vv)
    rr=rows('secondary/piston/metrics.csv')
    latex(5,[r'$N_\lambda$',r'$h/\lambda$',r'$E_2$',r'$E_\infty$'],[[r['cellsPerWavelength'],f"{float(r['h_over_lambda']):.4f}",f"{float(r['relL2']):.3e}",f"{float(r['relLinf']):.3e}"] for r in rr])
    latex(6,[r'$N_\lambda$',r'$h/\lambda$',r'$E_{2,|P|}^{ff}$'],[[r['cellsPerWavelength'],f"{float(r['h_over_lambda']):.4f}",f"{float(r['farField_pressureMagnitude_relL2']):.3e}"] for r in rr])
    rr=rows('secondary/gorkov/mesh-convergence.csv')
    latex(7,['Level',r'$h_s/a$','Segments','Faces','Force [fN]',r'Error [\%]'],[[r['level'],f"{float(r['h_over_a']):.5f}",r['segments'],r['n_surface_faces'],f"{float(r['numerical_force_y_N'])*1e15:.6f}",f"{float(r['difference_from_gorkov'])*100:.3f}"] for r in rr])
    report['reproduced_numerical_figures']=12;report['regenerated_tables']=7
    (out/'validation-report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
