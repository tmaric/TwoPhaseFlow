# Working notes for agents in this repository

## HPC job hygiene (Lichtenberg / HHLR, TU Darmstadt)

**Never cancel jobs account-wide.** `scancel -u tm83tomy`, `scancel --me` and
equivalents are forbidden: the account runs other work concurrently (leia
studies, Bayesian-optimisation orchestration, long queued jobs), and an
account-wide cancel destroys it with no warning and no record of what was lost.

Record the SLURM job IDs you submit and cancel only those. Snakemake's SLURM
executor logs every submission, so they are always recoverable:

```bash
grep -o 'SLURM jobid [0-9]*' run.log | awk '{print $NF}' | sort -u > my_jobs.txt
scancel $(cat my_jobs.txt)
```

Check what is running before acting: `squeue -u tm83tomy -o "%.10i %.9T %.30j"`.
Jobs you did not submit are not yours to touch.

**`pkill` needs the same care, for two reasons.** Use a pattern that cannot
match the wrapper running it — `pkill -f 'snakemak[e]'`, not `pkill -f
snakemake` — or an ssh command whose own command line contains the pattern
kills its own session.

More importantly, *a pattern kill is account-wide just as `scancel -u` is*.
The account runs other Snakemake studies (`leia-curvature`, the Bayesian
optimisation workflows) from the same login nodes, and `pkill -f 'snakemak[e]'`
would end those drivers too, leaving their submitted jobs orphaned in the queue
with nobody to collect the results. Kill by recorded PID, or scope the pattern
to this repository's directory:

```bash
pgrep -af 'snakemak[e]' | grep TwoPhaseFlow | awk '{print $1}' > my_pids.txt
kill $(cat my_pids.txt)
```

Before killing anything, check whose it is: `readlink /proc/<pid>/cwd`.

**Do not modify `/work/scratch/tm83tomy/leia`.** It is a reference for SLURM and
Snakemake conventions (see `leia/CLUSTER.md` and `leia/profiles/slurm/config.yaml`)
and is read-only for work in this repository.

**Keep builds separate.** Use a dedicated `WM_PROJECT_USER_DIR` so an OpenFOAM
build here never overwrites another project's binaries; this repository's
verification workflow uses `$HOME/OpenFOAM/tpf-rfv-v2512`.

## Verification workflow

`verification/reversed-vortex/` reproduces the reversed-vortex advection study
with one command (`./run.sh`). See its README for the protocol, and note that
OpenFOAM's `adjustableRunTime` write control snaps the time step onto the write
times, so the write interval is part of the protocol, not a passive observer.
