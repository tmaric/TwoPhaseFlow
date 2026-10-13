#!/bin/bash
# The sphere-translation ladder end to end: generate, run, collect (the shared
# runner ../runAll-common.sh: --jobs K, --cases FILE, --case DIR, --gen, the
# rest to genCases.py).
#     source $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc && . <repo>/etc/davof-env.sh
#     ./runAll.sh --N 32 64                     # the local ladder, one case at a time
#     ./runAll.sh --jobs 24                     # the cluster list (36 cases: 12 rows x 3 N)
#     ./runAll.sh --case Cases/davof-quadraticFaces-exact/N32
cd "${0%/*}" || exit 1
BENCH_GLOB='Cases/*/N*' exec bash ../runAll-common.sh "$@"
