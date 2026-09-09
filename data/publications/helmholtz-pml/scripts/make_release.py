#!/usr/bin/env python3
"""Build a checked, deterministic ZIP for transferring the secondary-data deposit."""
from pathlib import Path
import argparse, hashlib, json, sys, zipfile
sys.dont_write_bytecode=True
from validate import ROOT, run

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    output=args.output.resolve()
    if output.is_relative_to(ROOT):raise ValueError('The export must be outside the immutable payload.')
    if output.exists():raise FileExistsError(output)
    result=run()
    output.parent.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((ROOT/'checksums.json').read_text())
    names=sorted(list(manifest['files'])+['checksums.json'])
    # TUdatalib recommends uncompressed containers when preserving directories.
    with zipfile.ZipFile(output,'x',compression=zipfile.ZIP_STORED) as z:
        for name in names:
            info=zipfile.ZipInfo('helmholtz-pml/'+name,date_time=(2026,9,9,0,0,0))
            info.compress_type=zipfile.ZIP_STORED;info.external_attr=0o100644<<16
            z.writestr(info,(ROOT/name).read_bytes())
    print(json.dumps(dict(archive=str(output),files=len(names),bytes=output.stat().st_size,sha256=hashlib.sha256(output.read_bytes()).hexdigest(),validation=result['status']),indent=2))

if __name__=='__main__':main()
