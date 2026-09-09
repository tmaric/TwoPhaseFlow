#!/usr/bin/env bash
#SBATCH --account=special00004
#SBATCH --job-name=acoustic-interface-release
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=3000
#SBATCH --time=00:30:00
set -o pipefail
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv python/3.12.12-vzbp
source /home/tm83tomy/OpenFOAM/OpenFOAM-v2512/etc/bashrc
repo=/work/scratch/tm83tomy/TwoPhaseFlowAcoustics
cd "$repo"
source scripts/bashrc
work="$repo/.review-testing/interface-implementation"
old="$repo/.review-testing/openfoam-user/platforms/$WM_OPTIONS"
export LIBRARY_PATH="$old/lib:${LIBRARY_PATH:-}"
export FOAM_USER_APPBIN="$work/build/bin" FOAM_USER_LIBBIN="$work/build/lib"
export PATH="$FOAM_USER_APPBIN:$old/bin:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$old/lib:$LD_LIBRARY_PATH"
export CPATH="$(mpicc --showme:incdirs | tr ' ' ':')${CPATH:+:$CPATH}"
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
set -e
mkdir -p "$FOAM_USER_APPBIN" "$FOAM_USER_LIBBIN"
BUILD_JOBS=8 bash testsuite/acousticInterface/build.sh > "$work/release-build-$SLURM_JOB_ID.log" 2>&1
source "$work/venv/bin/activate"
python3 -m pytest --confcutdir=testsuite/acousticInterface testsuite/acousticInterface/test_invariants.py \
 --basetemp="$work/pytest-release-$SLURM_JOB_ID" -q > "$work/pytest-release-$SLURM_JOB_ID.log" 2>&1
python3 testsuite/acousticInterface/compare.py --output "$work/release-gate-$SLURM_JOB_ID" \
 --sizes 160 --offsets 0.25 --ranks 1,8 --baseline "$work/acousticHelmholtzFoam-baseline" \
 > "$work/release-gate-$SLURM_JOB_ID.log" 2>&1
python3 testsuite/acousticInterface/audit_parallel.py \
 --case "$work/study-54503635-2-sphere/N48-a0-transmission-np8" > "$work/parallel-audit-$SLURM_JOB_ID.log" 2>&1
printf 'RELEASE_CHECKS_COMPLETED\n'
