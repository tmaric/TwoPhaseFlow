#!/usr/bin/env python3
"""Create a compact, reproducible report from retained comparison outputs."""
import argparse,json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read(path): return json.loads(path.read_text())
def order(rows):
    rows=sorted(rows,key=lambda x:x["N"])[-3:]
    return float(-np.polyfit(np.log([x["N"] for x in rows]),
                            np.log([x["pressureRelL2"] for x in rows]),1)[0])
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--results",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    cluster=a.results/"lichtenberg"
    def study(suffix):
        paths=sorted(cluster.glob(f"*{suffix}/results.json"))
        return paths[-1] if paths else None
    forward=study("-forward") or a.results/"layered-forward/results.json"
    reverse=study("-reverse")
    sphere=study("-sphere")
    if not sphere or len(read(sphere))<24: sphere=a.results/"sphere-refinement-local/results.json"
    oblique=study("-oblique"); homogeneous=study("-homogeneous")
    wedges=sorted(cluster.glob("wedge-*/results.json"))
    wedge=wedges[-1] if wedges else None
    regular=sorted(cluster.glob("regular-oblique-*/results.json"))
    regular=regular[-1] if regular else a.results/"oblique-regular/results.json"
    axisymmetric=a.results/"sphere-wedge/results.json"
    sphereRows=read(sphere); forwardRows=read(forward)
    summary={"sphereSource":str(sphere),"forwardSource":str(forward),"studies":{}}
    for title,path in [("forward",forward),("reverse",reverse),("sphere",sphere),
                       ("oblique",oblique),("homogeneous",homogeneous),("wedge",wedge),
                       ("regularOblique",regular if regular.exists() else None),
                       ("axisymmetricSphere",axisymmetric if axisymmetric.exists() else None)]:
        if path:
            rows=read(path)
            summary["studies"][title]={
                "runs":len(rows),"failedRuns":sum(not x.get("completed",True) for x in rows),
                "rows":rows}
            supplement=path.parent/"metrics-supplement.json"
            if supplement.exists(): summary["studies"][title]["supplementalMetrics"]=read(supplement)
    labels={"legacy":"Legacy","geometry":"Geometric areas + legacy flux",
            "transmission":"Geometric areas + transmission flux"}
    colours={"legacy":"#333333","geometry":"#4387ad","transmission":"#b44926"}
    fig,axes=plt.subplots(1,2,figsize=(10,4.2),layout="constrained")
    for ax,rows,title in [(axes[0],forwardRows,"Layered interface: quarter-cell offset"),
                          (axes[1],sphereRows,"Penetrable sphere: 3D")]:
        for mode in labels:
            selected=sorted([x for x in rows if x["mode"]==mode and x["ranks"]==1
                             and x.get("offset",.25)==.25 and x.get("completed",True)],
                            key=lambda x:x["N"])
            ax.loglog([x["N"] for x in selected],[x["pressureRelL2"] for x in selected],
                      marker={"legacy":"o","geometry":"s","transmission":"^"}[mode],
                      ls="--" if mode=="geometry" else "-",markersize=4,
                      color=colours[mode],label=labels[mode])
        ax.set_title(title,fontsize=11);ax.set_xlabel("Cells along domain")
        ax.set_ylabel("Relative complex-pressure L2 error")
        ax.grid(True,which="both",alpha=.2)
    axes[1].legend(fontsize=8,loc="lower left")
    fig.savefig(a.output/"pressure-comparison.svg",metadata={"Date":None})
    fig.savefig(a.output/"pressure-comparison.png",dpi=180)
    plt.close(fig)
    summary["sphereOrders"]={m:order([x for x in sphereRows if x["mode"]==m and x["ranks"]==1]) for m in labels}
    summary["sphereCost"]={}
    for mode in labels:
        rows=[x for x in sphereRows if x["mode"]==mode and x["N"]==48]
        row=max(rows,key=lambda x:x["ranks"])
        log=sphere.parent/f"N48-a0-{mode}-np{row['ranks']}"/"log.solver"
        text=log.read_text()
        memory=re.findall(r"INFOG\(22\).*:\s*(\d+)",text)
        nonzeros=re.findall(r"total: nonzeros=(\d+)",text)
        summary["sphereCost"][mode]={"ranks":row["ranks"],"seconds":row["seconds"],
            "mumpsFactorMemorySumMB":int(memory[-1]) if memory else None,
            "matrixNonzeros":int(nonzeros[-1]) if nonzeros else None}
    (a.output/"verification-summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False))
    print(json.dumps({k:v for k,v in summary.items() if k!="studies"},indent=2))
if __name__=="__main__":main()
