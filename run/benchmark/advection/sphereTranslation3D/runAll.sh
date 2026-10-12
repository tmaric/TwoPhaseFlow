#!/bin/bash
# The sphere-translation ladder end to end: generate, run (coarse to fine, one
# case at a time, timed), collect.
#     source $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc && . <repo>/etc/davof-env.sh
#     ./runAll.sh --N 32 64              # the local ladder
#     ./runAll.sh                        # N = 32, 64, 128 (128: 4.2e6 cells)
cd "${0%/*}" || exit 1
python3 genCases.py "$@" || exit 1
for n in $(ls -d Cases/*/N* | sed 's:.*/N::' | sort -n | uniq); do
    for c in Cases/*/N$n/; do
        t0=$(date +%s)
        (cd "$c" && ./Allrun > log.Allrun 2>&1)
        rc=$?
        echo "$c rc=$rc  $(( $(date +%s) - t0 )) s"
    done
done
python3 getResults.py
