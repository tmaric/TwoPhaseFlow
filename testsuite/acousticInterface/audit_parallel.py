#!/usr/bin/env python3
"""Check processor-face operator cancellation from compiled solver diagnostics."""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from compare import field

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--case",type=Path,required=True)
    a=ap.parse_args()
    ranks=sorted(a.case.glob("processor[0-9]*"),key=lambda p:int(p.name[9:]))
    if len(ranks)<2: raise ValueError("A decomposed run with raw processor fields is required")
    counts=[len(field(p/"1/Pre")) for p in ranks]
    ends=np.cumsum(counts)
    owner_rank=lambda i:int(np.searchsorted(ends,i,side="right"))
    pairs=defaultdict(list)
    max_ranks=0
    for rank,path in enumerate(ranks):
        groups={}
        with (path/"postProcessing/acousticInterface/1/operators.tsv").open() as stream:
            for row in csv.DictReader(stream,delimiter="\t"):
                f=int(row["face"])
                if f not in groups:
                    groups[f]=(int(row["ownerGlobal"]),int(row["neighbourGlobal"]),{})
                groups[f][2][int(row["columnGlobal"])]=float(row["weight"])
        for o,n,weights in groups.values():
            max_ranks=max(max_ranks,len({owner_rank(k) for k in weights}))
            if owner_rank(o)!=owner_rank(n): pairs[tuple(sorted((o,n)))].append(weights)
    if not pairs: raise ValueError("No affected processor faces were exercised")
    worst=0
    for key,operators in pairs.items():
        if len(operators)!=2: raise ValueError(f"Missing/duplicate face operator: {key}")
        x,y=operators
        if x.keys()!=y.keys(): raise ValueError(f"Different processor stencils: {key}")
        scale=sum(abs(v) for v in x.values())
        defect=max(abs(x[k]+y[k]) for k in x)/max(scale,1e-300)
        worst=max(worst,defect)
    report={"processorFaces":len(pairs),"maxStencilRanks":max_ranks,
            "sharedOperatorDefect":worst,"passed":bool(worst<=1e-12)}
    (a.case/"parallel-operator-audit.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    if not report["passed"]: raise SystemExit(1)
if __name__=="__main__": main()
