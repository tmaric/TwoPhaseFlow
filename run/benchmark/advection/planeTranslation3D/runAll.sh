#!/bin/bash
# The plane-translation gate end to end: generate, run, collect (the shared
# runner ../runAll-common.sh: --jobs K, --cases FILE, --case DIR, --gen, the
# rest to genCases.py).
#     source $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc && . <repo>/etc/davof-env.sh
#     ./runAll.sh                               # the cluster list, one case at a time
#     ./runAll.sh --jobs 8 --N 20 40            # eight at a time
#     ./runAll.sh --case Cases/normal/davof-detrixheAslam-Youngs0.1/N20
cd "${0%/*}" || exit 1
BENCH_GLOB='Cases/*/*/N*' exec bash ../runAll-common.sh "$@"
