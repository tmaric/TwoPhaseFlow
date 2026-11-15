#!/usr/bin/env bash
#SBATCH --account=special00004
#SBATCH --job-name=cpc261115-build
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=4000
#SBATCH --time=00:30:00
#SBATCH --output=/work/scratch/tm83tomy/ultrasound-cpc-20261115/logs/build-%j.log
source /work/scratch/tm83tomy/ultrasound-cpc-20261115/scripts/hpc-env.sh
set -eo pipefail
cd "$CPC_REPO"
test "$(git rev-parse HEAD)" = 7f00fc55a67aefbcfa88bdb46a4fc378bac0db0c
mkdir -p "$FOAM_USER_APPBIN" "$FOAM_USER_LIBBIN"
wmake -j 8 libso src/VoF
BUILD_JOBS=8 bash testsuite/acousticInterface/build.sh
wmake -j 8 apps/utilities/mesh/setPMLFields
bash apps/utilities/preProcessing/setAlphaField/Allwmake
command -v acousticHelmholtzFoam setPMLFields setAlphaField cartesian2DMesh
test "$(command -v acousticHelmholtzFoam)" = "$FOAM_USER_APPBIN/acousticHelmholtzFoam"
ldd "$FOAM_USER_APPBIN/acousticHelmholtzFoam" > "$CPC_ROOT/logs/solver-dependencies.txt"
! grep -q 'not found' "$CPC_ROOT/logs/solver-dependencies.txt"
sha256sum "$FOAM_USER_APPBIN/acousticHelmholtzFoam" "$FOAM_USER_APPBIN/setPMLFields" "$FOAM_USER_LIBBIN/libacousticInterface.so" > "$CPC_ROOT/logs/build-sha256.txt"
env | grep -E '^(WM_PROJECT|FOAM_USER|PETSC|TPF_PROJECT|SLURM_JOB_ID)' > "$CPC_ROOT/logs/build-environment.txt"
python3 -m pytest --confcutdir=testsuite/acousticInterface testsuite/acousticInterface/test_invariants.py --basetemp="$CPC_ROOT/results/pytest-build" -q
printf '\nCPC_BUILD_COMPLETE\n'
