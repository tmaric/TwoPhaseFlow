"""Exact kinematic illustration. No ODE integration or CFD simulation.

Run: python3 two_motions_storyboard.py [OUTPUT_DIRECTORY] [--animate]
Outputs a vector PDF, SVG and PNG preview; --animate also creates a GIF.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch, Circle
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox

OUT=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent
OUT.mkdir(exist_ok=True,parents=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'svg.fonttype':'none'})
INK='#20323e'; MUTED='#52616a'; GRID='#e3e8eb'; PHASE='#377c7b'; FRAME='#586b7d'
MARKERS=['#c4480a','#1462b0']
MARKER_SHAPES=['o','D']
CENTER=np.array([.55,0.]);RADIUS=.26;SHEAR=1.6
angles=np.linspace(0,2*np.pi,361)
initial=CENTER+RADIUS*np.c_[np.cos(angles),np.sin(angles)]
labels=CENTER+RADIUS*np.c_[np.cos([np.pi/3,7*np.pi/6]),np.sin([np.pi/3,7*np.pi/6])]

def poses(s, points):
    # Force exact endpoint poses for drawing, without hiding a numerical solver error.
    k=0. if s in (0.,1.) else SHEAR*np.sin(np.pi*s)**2
    S=np.array([[1.,k],[0.,1.]])
    angle=2*np.pi*s
    Q=np.eye(2) if s in (0.,1.) else np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    base=CENTER+(points-CENTER)@S.T
    frame=points@Q.T
    combined=base@Q.T
    return base,frame,combined,Q

def arrow(ax,p,q,color,width=1.0,size=8,**kw):
    ax.add_patch(FancyArrowPatch(p,q,arrowstyle='-|>',mutation_scale=size,lw=width,color=color,**kw))

def scene(ax,s,row):
    ax.set(xlim=(-1.12,1.12),ylim=(-1.12,1.12),aspect='equal')
    ax.set_xticks([]);ax.set_yticks([])
    for sp in ax.spines.values():sp.set_color('#d7e0e5');sp.set_linewidth(.65)
    ax.axhline(0,color=GRID,lw=.65,zorder=0);ax.axvline(0,color=GRID,lw=.65,zorder=0)
    if row>0:
        # The orbit is a reference for the rigid frame, not the deforming boundary.
        ax.add_patch(Circle((0,0),CENTER[0],fill=False,color='#bdcee4',ls=(0,(2,2)),lw=.65,zorder=0))
    shape=poses(s,initial)[row]
    points=poses(s,labels)[row]
    ax.add_patch(Polygon(shape,closed=True,fc=('#f7f8f9' if row==1 else '#f4faf9'),
                         ec=('#6c7b86' if row==1 else PHASE),lw=1.05,zorder=2))
    for i,p in enumerate(points):
        ax.plot(*p,marker=MARKER_SHAPES[i],ms=5.7,mec='white',mew=.65,color=MARKERS[i],zorder=5)
    ax.plot(0,0,'+',ms=5,mew=.85,color=INK,zorder=6)
    if row>0:
        Q=poses(s,initial)[3]
        for vec in [np.array([.20,0.]),np.array([0.,.16])]:
            end=Q@vec
            arrow(ax,(0,0),end,FRAME,width=1.25,size=7,zorder=4)
    if row==1 and s==0.:
        ax.text(-.30,-.24,'$c$',fontsize=10,color=INK)
        th=np.linspace(.16,np.pi*.47,30)
        curve=.88*np.c_[np.cos(th),np.sin(th)]
        ax.plot(*curve.T,color=FRAME,lw=1.)
        arrow(ax,curve[-3],curve[-1],FRAME,size=7)

def storyboard():
    fig=plt.figure(figsize=(7.4,7.1),facecolor='white')
    gs=fig.add_gridspec(3,5,left=.035,right=.985,top=.79,bottom=.18,hspace=.50,wspace=.095)
    ts=[0.,.25,.5,.75,1.]
    times=['0',r'$T/4$',r'$T/2$',r'$3T/4$',r'$T$']
    headings=['A  Base deformation',
              'B  Frame motion alone',
              'C  Combined motion']
    for row in range(3):
        for col,(s,t) in enumerate(zip(ts,times)):
            ax=fig.add_subplot(gs[row,col]);scene(ax,s,row)
            if row==0:
                ax.text(.5,1.025,t,transform=ax.transAxes,ha='center',va='bottom',fontsize=10.5,color=INK)
            if row==2:
                label=['initial\ndisc','deformed\nand turned','largest\nbase shear','same base shape\nas at T/4','exact\ninitial state'][col]
                ax.text(.5,-.12,label,transform=ax.transAxes,ha='center',va='top',fontsize=9,color=MUTED)
        box=gs[row,0].get_position(fig)
        offset=.054 if row==0 else .024
        fig.text(.035,box.y1+offset,headings[row],fontsize=10,color=INK)
    handles=[Line2D([],[],color=MARKERS[0],marker=MARKER_SHAPES[0],lw=0,ms=6,label='orange material point'),
             Line2D([],[],color=MARKERS[1],marker=MARKER_SHAPES[1],lw=0,ms=6,label='blue material point'),
             Line2D([],[],color=FRAME,lw=1.5,label='carried frame axes'),
             Line2D([],[],color=GRID,lw=1.5,label='fixed laboratory axes')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.51,.012),ncol=2,frameon=False,fontsize=9.5,handlelength=1.4,columnspacing=2.1)
    base=OUT/'two_motions_storyboard'
    for ext in ('pdf','svg','png'):fig.savefig(str(base)+'.'+ext,dpi=210,bbox_inches=Bbox.from_extents(0,0,7.4,6.22))
    plt.close(fig)

def animation():
    from PIL import Image
    fig,axs=plt.subplots(1,3,figsize=(9.0,3.4))
    fig.subplots_adjust(left=.025,right=.985,top=.78,bottom=.08,wspace=.15)
    frames=[]
    times=np.r_[np.zeros(5),np.linspace(0.,1.,81),np.ones(10)]
    for s in times:
        for row,ax in enumerate(axs):
            ax.clear();scene(ax,float(s),row)
            ax.set_title(['A  Base deformation','B  Frame motion alone','C  Both motions together'][row],fontsize=10.5,color=INK,pad=7)
        fig.suptitle(f'Two returning motions  |  t / T = {s:.2f}',fontsize=14,fontweight='bold',color=INK,y=.96)
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3]).convert('RGB'))
    frames[0].save(OUT/'two_motions_animation.gif',save_all=True,append_images=frames[1:],duration=75,loop=0,optimize=False)
    plt.close(fig)

if __name__=='__main__':
    assert np.max(np.abs(poses(0.,initial)[2]-initial))==0.
    assert np.max(np.abs(poses(1.,initial)[2]-initial))==0.
    assert np.allclose(poses(.25,initial)[0],poses(.75,initial)[0])
    assert np.allclose(poses(.25,initial)[2],-poses(.75,initial)[2])
    storyboard()
    if '--animate' in sys.argv:animation()
    print(OUT/'two_motions_storyboard.pdf')
