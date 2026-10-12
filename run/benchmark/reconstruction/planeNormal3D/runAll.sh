#!/bin/bash
# The reconstruction ladder of this benchmark end to end: generate, run, collect.
#     source $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc && . <repo>/etc/davof-env.sh
#     ./runAll.sh [genCases.py arguments]
# Cases run coarse to fine (an alphabetical glob would start with N160, the
# 4.1e6-cell cutCellImpFunc cases needing ~8 GB each), one at a time, timed.
cd "${0%/*}" || exit 1
python3 genCases.py "$@" || exit 1
for n in $(ls -d Cases/*/*/N* | sed 's:.*/N::' | sort -n | uniq); do
    for c in Cases/*/*/N$n/; do
        t0=$(date +%s)
        (cd "$c" && ./Allrun > log.Allrun 2>&1)
        rc=$?
        echo "$c rc=$rc  $(( $(date +%s) - t0 )) s"
    done
done
python3 getResults.py
