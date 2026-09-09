#!/usr/bin/env bash
#SBATCH --account=special00004
#SBATCH --job-name=acoustic-interface-wedge
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
export FOAM_MODULE_APPBIN="$FOAM_USER_APPBIN" FOAM_MODULE_LIBBIN="$FOAM_USER_LIBBIN"
export LIBRARY_PATH="$old/lib"
export PATH="$FOAM_USER_APPBIN:$old/bin:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$old/lib:$LD_LIBRARY_PATH"
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PETSC_OPTIONS="-ksp_view"
source "$work/venv/bin/activate"
set -e
# Build the source snapshot used for the WSL cfMesh tests. Some cluster
# OpenFOAM installations ship an empty plugin placeholder, not the sources.
cfmeshSource="$work/cfmesh-source"
if ! test -f "$cfmeshSource/meshLibrary/Make/files"; then
    echo "Supply the cfMesh source snapshot in $cfmeshSource before submission" >&2
    exit 2
fi
wmake -j8 libso "$cfmeshSource/meshLibrary" > "$work/cfmesh-build.log" 2>&1
wmake -j8 "$cfmeshSource/executables/cartesian2DMesh" >> "$work/cfmesh-build.log" 2>&1
python3 testsuite/acousticInterface/wedge.py --output "$work/wedge-$SLURM_JOB_ID" \
    --sizes 0.0004,0.0002 --aspects 1,2 --ranks 1,8
printf 'INTERFACE_WEDGE_STUDY_COMPLETED\n'
