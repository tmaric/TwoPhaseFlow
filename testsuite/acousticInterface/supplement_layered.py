#!/usr/bin/env python3
"""Add reflection/transmission phases and velocity errors from retained fields."""
import argparse,json
from pathlib import Path
from compare import layered_metrics

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--study",type=Path,required=True)
    ap.add_argument("--reverse",action="store_true")
    a=ap.parse_args()
    rows=json.loads((a.study/"results.json").read_text())
    supplemental=[]
    for row in rows:
        case=a.study/f"N{row['N']}-o{row['offset']}-{row['mode']}-np{row['ranks']}"
        _,metrics=layered_metrics(case,row["N"],row["offset"],a.reverse)
        if abs(metrics["pressureRelL2"]-row["pressureRelL2"])>1e-12:
            raise RuntimeError(f"Retained pressure results changed: {case}")
        supplemental.append({k:row[k] for k in ["N","offset","mode","ranks"]}|metrics)
    dest=a.study/"metrics-supplement.json"
    if dest.exists(): raise FileExistsError(dest)
    dest.write_text(json.dumps(supplemental,indent=2,allow_nan=False))
    print(f"{len(supplemental)} retained cases audited: {dest}")
if __name__=="__main__": main()
