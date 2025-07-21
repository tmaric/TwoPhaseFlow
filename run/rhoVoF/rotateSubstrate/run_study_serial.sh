#!/usr/bin/env bash

PARAM_FILE_PREFIX=${1%%.parameter}
SOLVER=$2
TEMPLATECASE="${3:-templateCase}"
# N_CASES_IN_PARALLEL=$3

create-parameter-study.py -p "$SOLVER" -t "$TEMPLATECASE" "$PARAM_FILE_PREFIX".parameter  && \
#bulkeval "$SOLVER"-"$PARAM_FILE_PREFIX"_000 "sed -i '1 i\#!/bin/bash' isoAdv.sbatch"
# bulkeval "$SOLVER"-"$PARAM_FILE_PREFIX"_000 "sed -i '1 i\#!/bin/bash' hexrefinedMesh.sh"
# bulkeval "$SOLVER"-"$PARAM_FILE_PREFIX"_000 "sed -i '1 i\#!/bin/bash' initFields.sh"
initialize-parameter-study.py "$SOLVER"-"$PARAM_FILE_PREFIX"_000 -m blockMesh -f initFields.sh -j

