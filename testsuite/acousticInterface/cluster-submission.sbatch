#!/usr/bin/env bash
#SBATCH --account=special00004
#SBATCH --job-name=acoustic-submission
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=4000
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
export FOAM_USER_APPBIN="$work/build/bin" FOAM_USER_LIBBIN="$work/build/lib"
export PATH="$FOAM_USER_APPBIN:$old/bin:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$old/lib:$LD_LIBRARY_PATH"
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PETSC_OPTIONS="-ksp_view"
source "$work/venv/bin/activate"
set -e
out="$repo/.review-testing/submission-geometric-averaging"
mkdir -p "$out"
python3 -m pytest --confcutdir=testsuite/acousticInterface testsuite/acousticInterface/test_invariants.py \
 --basetemp="$out/pytest-$SLURM_JOB_ID" -q > "$out/pytest-$SLURM_JOB_ID.log" 2>&1
python3 testsuite/acousticInterface/submission_check.py --output "$out/pml-$SLURM_JOB_ID" \
 --baseline "$work/acousticHelmholtzFoam-baseline" --ranks 1,8 > "$out/pml-$SLURM_JOB_ID.log" 2>&1
printf 'SUBMISSION_CHECKS_COMPLETED\n'
