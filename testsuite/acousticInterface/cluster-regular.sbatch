#!/usr/bin/env bash
#SBATCH --account=special00004
#SBATCH --job-name=acoustic-interface-regular
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
python3 testsuite/acousticInterface/advanced.py --output "$work/regular-oblique-$SLURM_JOB_ID" \
 --scenario oblique --mesh-style regular --sizes 48,96,192 --angles 15,45,75 --ranks 1,8
