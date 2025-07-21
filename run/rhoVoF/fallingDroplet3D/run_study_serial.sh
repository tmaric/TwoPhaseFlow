#!/usr/bin/env bash

PARAM_FILE_PREFIX=${1%%.parameter}
SOLVER=$2
TEMPLATECASE="${3:-templateCase}"
#TEMPLATE=$3

create-parameter-study.py  -t "$TEMPLATECASE" -p "$SOLVER" "$PARAM_FILE_PREFIX".parameter && \
initialize-parameter-study.py "$SOLVER"-"$PARAM_FILE_PREFIX"_000 -m blockMesh -f initFields.sh 
