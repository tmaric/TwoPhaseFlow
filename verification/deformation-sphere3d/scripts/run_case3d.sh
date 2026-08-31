#!/bin/bash
# Mesh, initialise and advect one materialised case.
#
#   <script> <case-dir> <np> <reconstruct> [mpi launcher ...]
#
# RECONSTRUCT is 0 for the convergence cases and 1 for the snapshot cases.
# volumeFractionError is a function object and writes correct global values from
# a parallel run without reconstruction, so only the snapshots -- which need the
# volume fraction itself, on the whole mesh, at several times -- pay for
# reconstructPar. At 2048^2 and 256^3 that is the difference between minutes and
# hours.
#
# COMPLETION IS JUDGED FROM THE SOLVER LOG, NOT THE EXIT CODE. This OpenFOAM
# build intermittently corrupts the heap while tearing down static objects,
# after the time loop has finished and everything has been written:
#
#     End
#     malloc_consolidate(): invalid chunk size
#
# Observed on both serial and parallel runs, 2D and 3D. Trusting the exit code
# discards complete, correct results, so instead we require that the solver
# printed "End" and that the error functional wrote its file.
set -e
CASE="$1"; NP="$2"; RECON="$3"; shift 3; LAUNCHER="$*"
cd "$CASE"

blockMesh                 > log.blockMesh      2>&1
# restore0Dir is a RunFunctions shell *function*, not an executable
rm -rf 0 && cp -r 0.orig 0
setAlphaField             > log.setAlphaField  2>&1

set +e
if [ "$NP" -gt 1 ]; then
    decomposePar -force   > log.decomposePar   2>&1
    $LAUNCHER advectorVoF -parallel > log.advectorVoF 2>&1
    rc=$?
    if [ "$RECON" -eq 1 ]; then
        reconstructPar    > log.reconstructPar 2>&1
    fi
else
    advectorVoF           > log.advectorVoF    2>&1
    rc=$?
fi
set -e

if ! grep -q '^End$' log.advectorVoF; then
    echo "advectorVoF did not reach the end of the time loop (exit $rc)" >&2
    tail -20 log.advectorVoF >&2
    exit 1
fi
if [ "$rc" -ne 0 ]; then
    echo "NOTE: solver reached End but exited $rc -- known teardown crash in" >&2
    echo "      this build; the time loop and all output completed." >&2
fi

SRC=postProcessing/volumeFractionError/0/volumeFractionError.dat
[ -f "$SRC" ] || { echo "missing $SRC" >&2; exit 1; }
cp "$SRC" volumeFractionError.dat
