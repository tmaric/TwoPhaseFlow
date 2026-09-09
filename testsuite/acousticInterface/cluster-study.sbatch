#!/usr/bin/env bash
#SBATCH --account=special00004
#SBATCH --job-name=acoustic-interface-study
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=4000
#SBATCH --time=00:30:00
#SBATCH --array=0-2
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
out="$work/study-$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID"
status=0
case "$SLURM_ARRAY_TASK_ID" in
0)
 python3 testsuite/acousticInterface/compare.py --output "$out-forward" --ranks 1,8
 python3 testsuite/acousticInterface/compare.py --output "$out-reverse" --ranks 1,8 --reverse
 ;;
1)
 python3 testsuite/acousticInterface/advanced.py --output "$out-oblique" \
    --scenario oblique --sizes 48,96,192 --angles 15,45,75 --ranks 1,8 || status=1
 python3 testsuite/acousticInterface/advanced.py --output "$out-homogeneous" \
    --scenario homogeneous --sizes 48,96,192 --ranks 1,8
 ;;
2)
 python3 testsuite/acousticInterface/advanced.py --output "$out-sphere" \
    --scenario sphere --sizes 16,24,32,48 --ranks 1,8
 ;;
esac
printf 'INTERFACE_STUDY_COMPLETED status=%s\n' "$status"
exit "$status"
