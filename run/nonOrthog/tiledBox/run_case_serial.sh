#!/usr/bin/env bash

PARAM_FILE_PREFIX=${1%%.parameter}
SOLVER=$2
# N_CASES_IN_PARALLEL=$3

create-parameter-study.py -p "$SOLVER" -t template "$PARAM_FILE_PREFIX".parameter  && \
bulkeval "$SOLVER"-"$PARAM_FILE_PREFIX"_000 "sed -i '1 i\#!/bin/bash' interFlow.sbatch"
bulkeval "$SOLVER"-"$PARAM_FILE_PREFIX"_000 "sed -i '1 i\#!/bin/bash' initFields.sh"
initilize-parameter-study.py "$SOLVER"-"$PARAM_FILE_PREFIX"_000 -m blockMesh -f initFields.sh

