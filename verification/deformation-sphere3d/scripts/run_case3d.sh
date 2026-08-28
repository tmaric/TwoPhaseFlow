#!/bin/bash
# Mesh, initialise and advect one materialised 3D deformation case.
#   run_case3d.sh <case-dir> <np> [mpi launcher ...]
#
# No flux generator: movingFrameFlow3D builds phi from the analytic vector
# potential, which is what replaces generateUDeform + CorrectPhi.
set -e
CASE="$1"; NP="$2"; shift 2; LAUNCHER="$*"
cd "$CASE"
blockMesh                 > log.blockMesh      2>&1
rm -rf 0 && cp -r 0.orig 0
setAlphaField             > log.setAlphaField  2>&1
if [ "$NP" -gt 1 ]; then
    decomposePar -force   > log.decomposePar   2>&1
    $LAUNCHER advectorVoF -parallel > log.advectorVoF 2>&1
    reconstructPar        > log.reconstructPar 2>&1 || true
else
    advectorVoF           > log.advectorVoF    2>&1
fi
cp postProcessing/volumeFractionError/0/volumeFractionError.dat volumeFractionError.dat
