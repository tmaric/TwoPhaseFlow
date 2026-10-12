#!/bin/bash
# The plane-translation gate end to end: generate, run (one case at a time), collect.
#     source $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc && . <repo>/etc/davof-env.sh
#     ./runAll.sh [genCases.py arguments]
cd "${0%/*}" || exit 1
python3 genCases.py "$@" || exit 1
for c in Cases/*/*/N*/; do
    t0=$(date +%s)
    (cd "$c" && ./Allrun > log.Allrun 2>&1)
    rc=$?
    echo "$c rc=$rc  $(( $(date +%s) - t0 )) s"
done
python3 getResults.py
