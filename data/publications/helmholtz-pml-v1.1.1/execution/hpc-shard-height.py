from pathlib import Path
import csv,json,subprocess,math,shutil
root=Path('/work/scratch/tm83tomy/ultrasound-cpc-20261115')
assert not (root/'height-shard-plan.json').exists(),'Already planned'
plan=[]
for series in ['empty','1','2','3']:
 task=root/'results'/f'height-{series}'
 src=task/'case/sweepResults/heightSweepResults.csv'
 rows=list(csv.DictReader(src.open()))
 assert rows and len(rows)<100
 key='empty' if series=='empty' else 'ab'+series
 for i,r in enumerate(rows):
  assert int(r['index'])==i and r['series_key']==key
  assert abs(float(r['height_fac'])-(.46+.18*i/99))<2e-15
  assert math.isfinite(float(r['Fn_reflector1_axisym_N']))
 shutil.copy2(src,task/'seed.csv')
 (task/'SUPERSEDED').write_text('Original job intentionally stopped after completed rows were saved; remaining heights sharded. See height-shard-plan.json.\n')
 remaining=list(range(len(rows),100))
 width=math.ceil(len(remaining)/4)
 for part,start in enumerate(range(len(rows),100,width)):
  end=min(start+width,100)-1
  plan.append(dict(series=series,part=part,start=start,end=end,count=end-start+1,
     fac_min=format(.46+.18*start/99,'.16g'),fac_max=format(.46+.18*end/99,'.16g')))
(root/'height-shard-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
print(json.dumps(plan,indent=2))
