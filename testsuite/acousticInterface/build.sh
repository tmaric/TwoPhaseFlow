#!/usr/bin/env bash
# Source the chosen OpenFOAM environment first, then scripts/bashrc.
set -euo pipefail
: "${TPF_PROJECT_DIR:?Source scripts/bashrc first}"
cd "$TPF_PROJECT_DIR"
wmake -j "${BUILD_JOBS:-4}" libso src/acousticInterface
wmake -j "${BUILD_JOBS:-4}" solver/acousticHelmholtzFoam
wmake -j "${BUILD_JOBS:-4}" apps/benchmark/testAcousticInterface
