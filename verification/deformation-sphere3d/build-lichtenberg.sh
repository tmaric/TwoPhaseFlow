#!/bin/bash
# Build TwoPhaseFlow on Lichtenberg against $HOME/OpenFOAM/OpenFOAM-v2512.
#
# Uses its OWN WM_PROJECT_USER_DIR ($HOME/OpenFOAM/tpf-rfv-v2512) so it never
# mixes with leia's or with any other TwoPhaseFlow build on the machine.
# Idempotent: wmake skips what is already current, so re-running is cheap.
set -o pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"

module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv || exit 1
export WM_COMPILER_TYPE=system
export WM_MPLIB=SYSTEMOPENMPI
source "$HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc" || exit 1

export WM_PROJECT_USER_DIR="$HOME/OpenFOAM/tpf-rfv-v2512"
export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
export PATH="$FOAM_USER_APPBIN:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH"
export WM_NCOMPPROCS="${WM_NCOMPPROCS:-16}"

echo "repo     : $REPO  ($(cd "$REPO" && git log --oneline -1))"
echo "OF       : $WM_PROJECT_VERSION  ($WM_OPTIONS)"
echo "user dir : $WM_PROJECT_USER_DIR"

cd "$REPO" || exit 1
./Allwmake
rc=$?
echo "ALLWMAKE_EXIT=$rc"

# the case-local flux generator used by every case
cd "$REPO/testsuite/advection/test-vortexShearedDisc" && wmake generateUVortex2D
echo "generateUVortex2D: $(command -v generateUVortex2D || echo MISSING)"
exit $rc
