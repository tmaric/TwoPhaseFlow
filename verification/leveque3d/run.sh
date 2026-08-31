#!/bin/bash
# ---------------------------------------------------------------------------
# One command to reproduce the 3D LeVeque deformation verification study.
#
#     ./run.sh                # auto: SLURM if sbatch exists, else local
#     ./run.sh local          # force the local profile
#     ./run.sh slurm          # force the SLURM profile
#
# Everything downstream -- building the library, materialising the cases,
# submitting one job per case, aggregating the convergence table and rendering
# the interface figures -- is a Snakemake dependency of `rule all`.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"

PROFILE="${1:-auto}"
if [ "$PROFILE" = "auto" ]; then
    if command -v sbatch >/dev/null 2>&1; then PROFILE=slurm; else PROFILE=local; fi
fi
echo "profile: $PROFILE"

if [ "${SKIP_BUILD:-0}" != "1" ]; then
    if [ "$PROFILE" = "slurm" ]; then
        echo "=== building TwoPhaseFlow (Lichtenberg toolchain) ==="
        ./build-lichtenberg.sh
    else
        echo "=== skipping build; source OpenFOAM and run ./Allwmake yourself for a local run ==="
    fi
fi

echo "=== snakemake ==="
exec snakemake --workflow-profile "profiles/$PROFILE" "${@:2}"
