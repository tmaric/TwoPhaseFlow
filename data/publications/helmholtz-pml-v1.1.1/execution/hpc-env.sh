#!/usr/bin/env bash
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv python/3.12.12-vzbp
source /home/tm83tomy/OpenFOAM/OpenFOAM-v2606/etc/bashrc
export CPC_ROOT=/work/scratch/tm83tomy/ultrasound-cpc-20261115
export CPC_REPO="$CPC_ROOT/TwoPhaseFlow"
export WM_PROJECT_USER_DIR="$CPC_ROOT/openfoam-user"
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"
export PATH="$FOAM_USER_APPBIN:/home/tm83tomy/OpenFOAM/cfmesh-v2606/platforms/$WM_OPTIONS/bin:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:/home/tm83tomy/OpenFOAM/cfmesh-v2606/platforms/$WM_OPTIONS/lib:$LD_LIBRARY_PATH"
source "$CPC_REPO/scripts/bashrc"
export CPATH="$(mpicc --showme:incdirs | tr ' ' ':')${CPATH:+:$CPATH}"
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
source "$CPC_ROOT/venv/bin/activate"
