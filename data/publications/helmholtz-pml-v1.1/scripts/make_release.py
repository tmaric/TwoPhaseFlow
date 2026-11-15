#!/usr/bin/env python3
"""Validate then create a deterministic directory-preserving ZIP."""
from pathlib import Path
import argparse,subprocess,sys,zipfile
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(); p.add_argument('--output',required=True,type=Path); a=p.parse_args()
target=a.output.resolve()
if target.exists() or target.is_relative_to(root): raise SystemExit('Use a new ZIP path outside the package')
subprocess.run([sys.executable,str(root/'scripts/validate.py')],check=True)
target.parent.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_STORED) as z:
 for source in sorted(root.rglob('*')):
  if source.is_file() and '__pycache__' not in source.parts:
   info=zipfile.ZipInfo('helmholtz-pml-v1.1/'+source.relative_to(root).as_posix(),date_time=(2026,11,15,0,0,0))
   info.external_attr=0o100644<<16
   z.writestr(info,source.read_bytes())
print(target)
