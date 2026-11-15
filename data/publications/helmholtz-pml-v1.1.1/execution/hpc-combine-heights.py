from pathlib import Path
import csv,json,math,hashlib
root=Path('/work/scratch/tm83tomy/ultrasound-cpc-20261115')
plan=json.loads((root/'height-shard-plan.json').read_text())
pending=[p for p in plan if not (root/f'results/height-shards/{p["series"]}-{p["part"]}/COMPLETED').exists()]
if pending: raise SystemExit('Pending shards: '+str([(p['series'],p['part']) for p in pending]))
combined=[]
for series in ['empty','1','2','3']:
 task=root/f'results/height-{series}'
 rows=list(csv.DictReader((task/'seed.csv').open()))
 sources=[{'path':str(task/'seed.csv'),'rows':len(rows),'range':[0,len(rows)-1]}]
 reference_hashes=tuple(line.split()[0] for line in (task/'binary-sha256.txt').read_text().splitlines())
 for part in [p for p in plan if p['series']==series]:
  shard=root/f'results/height-shards/{series}-{part["part"]}'
  assert tuple(line.split()[0] for line in (shard/'binary-sha256.txt').read_text().splitlines())==reference_hashes
  path=shard/'case/sweepResults/heightSweepResults.csv'
  data=list(csv.DictReader(path.open()))
  assert len(data)==part['count'] and len(rows)==part['start']
  for j,r in enumerate(data):
   i=part['start']+j
   assert abs(float(r['height_fac'])-(.46+.18*i/99))<3e-15
   assert math.isfinite(float(r['Fn_reflector1_axisym_N']))
   r['index']=str(i); rows.append(r)
  sources.append({'path':str(path),'rows':len(data),'range':[part['start'],part['end']], 'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 assert len(rows)==100
 target=task/'case/sweepResults/heightSweepResults.csv'
 with target.open('w',newline='') as f:
  w=csv.DictWriter(f,list(rows[0])); w.writeheader(); w.writerows(rows)
 (task/'height-combination-provenance.json').write_text(json.dumps(sources,indent=2)+'\n')
 (task/'COMPLETED').write_text(f'CPC_CAMPAIGN_COMPLETE height-{series}; seed plus verified shards\n')
 combined.extend(rows)
out=root/'heightSweepResults-combined.csv'
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,list(combined[0])); w.writeheader(); w.writerows(combined)
print('Verified and combined all 400 points:',out)
