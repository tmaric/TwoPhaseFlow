# Toolchain for the Lichtenberg cluster (TU Darmstadt). Sourced, not executed.
#
# This lives in a file rather than in the profile's `env_preamble` string on
# purpose: snakemake expands $VARIABLES in config values against the *submitting*
# environment, so a preamble that defines a variable and then uses it in the same
# string gets the later reference expanded to empty. Here the expansion happens
# in the job's own shell, in order.
#
# Modules and MPI transport settings follow leia/CLUSTER.md.

module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_mtl=^ofi,psm2

export WM_COMPILER_TYPE=system
export WM_MPLIB=SYSTEMOPENMPI
source "$HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc"

# Own user directory, so this build never mixes with leia's or with any other
# TwoPhaseFlow build on the machine.
export WM_PROJECT_USER_DIR="$HOME/OpenFOAM/tpf-rfv-v2512"
export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
export PATH="$FOAM_USER_APPBIN:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH"
