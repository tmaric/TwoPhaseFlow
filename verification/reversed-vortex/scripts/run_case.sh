#!/bin/bash
# Mesh, initialise and advect one already-materialised reversed-vortex case.
#
#   run_case.sh <case-dir> <np> [mpi launcher ...]
#
# Kept out of the Snakefile so the convergence cases and the snapshot cases run
# through exactly the same code path.
set -e

CASE="$1"
NP="$2"
shift 2
LAUNCHER="$*"

cd "$CASE"

blockMesh                 > log.blockMesh      2>&1
# restore0Dir is a RunFunctions shell *function*, not an executable
rm -rf 0 && cp -r 0.orig 0
setAlphaField             > log.setAlphaField  2>&1
generateUVortex2D         > log.generateU      2>&1

if [ "$NP" -gt 1 ]; then
    decomposePar -force   > log.decomposePar   2>&1
    $LAUNCHER advectorVoF -parallel > log.advectorVoF 2>&1
    reconstructPar -latestTime > log.reconstructPar 2>&1 || true
else
    advectorVoF           > log.advectorVoF    2>&1
fi

cp postProcessing/volumeFractionError/0/volumeFractionError.dat volumeFractionError.dat
