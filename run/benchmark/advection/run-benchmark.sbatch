#!/bin/bash
# One advection benchmark (planeTranslation3D or sphereTranslation3D) on one
# Lichtenberg node: the cases of the benchmark's genCases.py list run
# JOBS at a time through runAll.sh --jobs, then getResults.py collects.
#
#     cd /work/scratch/tm83tomy/davof/TwoPhaseFlow/run/benchmark/advection
#     sbatch --export=ALL,REPO=/work/scratch/tm83tomy/davof/TwoPhaseFlow,BENCH=$PWD/sphereTranslation3D run-benchmark.sbatch
#     sbatch --export=ALL,REPO=...,BENCH=$PWD/planeTranslation3D,GENARGS="--N 20 40" run-benchmark.sbatch
#
# Variables (--export): REPO the TwoPhaseFlow clone built on the cluster
# (etc/davof-env.sh derives everything from its own location), BENCH the
# benchmark directory, JOBS the cases run at a time (default 24 = ntasks),
# GENARGS the genCases.py arguments (default: the benchmark's cluster list),
# CASES a case-list file for runAll.sh --cases (then genCases.py is skipped).
# The cases are serial; the sphere at N = 128 holds 4.2e6 cells (a few GB
# each, eight of them in the cluster list), hence 24 x 3600 MB (Lichtenberg's job_submit plugin refuses --mem; --mem-per-cpu is the only form).
#SBATCH -J davof-advection
#SBATCH -A special00004
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3600
#SBATCH -t 04:00:00
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
: "${REPO:?set REPO=<TwoPhaseFlow clone> with --export}"
: "${BENCH:?set BENCH=<benchmark directory> with --export}"
JOBS="${JOBS:-${SLURM_NTASKS:-24}}"

# OpenFOAM's etc/bashrc reads unset variables, so set -u only after the environment is sourced.
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
source "$REPO/etc/davof-env.sh" || exit 1
set -u

echo "run-benchmark: $(date) on $(hostname), REPO=$REPO, BENCH=$BENCH, JOBS=$JOBS"
echo "run-benchmark: $(cd "$REPO" && git rev-parse --short HEAD) $(cd "$REPO" && git branch --show-current)"
cd "$BENCH" || exit 1
if [ -n "${CASES:-}" ]; then
    ./runAll.sh --jobs "$JOBS" --cases "$CASES"
else
    # shellcheck disable=SC2086
    ./runAll.sh --jobs "$JOBS" ${GENARGS:-}
fi
rc=$?
python3 getResults.py
echo "run-benchmark: done $(date), runAll rc=$rc"
exit $rc
