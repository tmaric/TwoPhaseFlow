#!/usr/bin/env bash

PARAM_FILE_PREFIX=${1%%.parameter}
SOLVER=$2
# N_CASES_IN_PARALLEL=$3

create-parameter-study.py -p "$SOLVER" "$PARAM_FILE_PREFIX".parameter # && \
initialize-parameter-study.py "$SOLVER"-"$PARAM_FILE_PREFIX"_000 -m blockMesh -f initFields.sh # && \
# argo-run-study.py "$SOLVER" -d "$SOLVER"-"$PARAM_FILE_PREFIX"_000 -n "$N_CASES_IN_PARALLEL" -j
