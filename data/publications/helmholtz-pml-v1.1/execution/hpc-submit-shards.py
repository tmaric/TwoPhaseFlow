from pathlib import Path
import subprocess,json,re,csv
root=Path('/work/scratch/tm83tomy/ultrasound-cpc-20261115')
assert not (root/'height-shard-jobs.tsv').exists()
plan=json.loads((root/'height-shard-plan.json').read_text())
with (root/'height-shard-jobs.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t'); w.writerow(['series','part','job_id','start','end'])
 for p in plan:
  args=['sbatch','--parsable','--account=special00004',f'--job-name=cpc-height-{p["series"]}-{p["part"]}',
   '--nodes=1','--ntasks=8','--mem-per-cpu=6000','--time=01:00:00',
   f'--output={root}/logs/height-shard-{p["series"]}-{p["part"]}-%j.log',str(root/'scripts/hpc-height-shard.sh'),
   p['series'],str(p['part']),str(p['count']),p['fac_min'],p['fac_max']]
  out=subprocess.check_output(args,text=True)
  job=re.findall(r'^([0-9]+)(?:;[^\s]+)?$',out,re.M)[-1]
  w.writerow([p['series'],p['part'],job,p['start'],p['end']]); f.flush()
  print(p['series'],p['part'],job,flush=True)
