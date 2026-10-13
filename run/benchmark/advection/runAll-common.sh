#!/bin/bash
# The shared runner of the advection benchmarks, called by <benchmark>/runAll.sh
# from the benchmark directory with BENCH_GLOB set to the case pattern
# (planeTranslation3D: 'Cases/*/*/N*', sphereTranslation3D: 'Cases/*/N*').
#
#     source $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc && . <repo>/etc/davof-env.sh
#     ./runAll.sh [--jobs K] [--cases FILE] [--case DIR ...] [--gen] [genCases.py arguments]
#
# Without a case list: genCases.py <arguments> generates the cases and every
# case matching BENCH_GLOB is run. With --cases FILE (one case directory per
# line, relative to the benchmark directory, '#' comments) or --case DIR the
# listed cases are run and genCases.py is skipped unless --gen is given (then
# it runs with the remaining arguments first). --jobs K runs K cases at a
# time (xargs -P; the cases are serial), the finest first when K > 1, the
# coarsest first when K = 1. Every case runs ./Allrun > log.Allrun; one line
# per case reports the exit code and the wall time. getResults.py collects.
set -u
[ -n "${BENCH_GLOB:-}" ] || { echo "runAll-common.sh: BENCH_GLOB is not set" >&2; exit 1; }

jobs=1
gen=0
cases=()
genargs=()
while [ $# -gt 0 ]; do
    case "$1" in
        --jobs)   jobs="$2"; shift 2 ;;
        --jobs=*) jobs="${1#--jobs=}"; shift ;;
        --gen)    gen=1; shift ;;
        --case)   cases+=("$2"); shift 2 ;;
        --cases)
            while IFS= read -r line; do
                line="${line%%#*}"; line="${line// /}"
                [ -n "$line" ] && cases+=("$line")
            done < "$2"
            shift 2 ;;
        *) genargs+=("$1"); shift ;;
    esac
done

if [ ${#cases[@]} -eq 0 ] || [ $gen -eq 1 ]; then
    python3 genCases.py "${genargs[@]}" || exit 1
fi
if [ ${#cases[@]} -eq 0 ]; then
    for c in $BENCH_GLOB; do
        [ -d "$c" ] && cases+=("$c")
    done
fi
[ ${#cases[@]} -gt 0 ] || { echo "runAll: no cases"; exit 1; }

# Order by N: the finest first for a parallel run, the coarsest first serially.
order="-n"
[ "$jobs" -gt 1 ] && order="-rn"
mapfile -t cases < <(
    for c in "${cases[@]}"; do
        c="${c%/}"
        echo "${c##*/N} $c"
    done | sort $order -k1,1 -s | cut -d' ' -f2-
)

run_one() {
    local c="${1%/}" t0 rc
    [ -x "$c/Allrun" ] || { echo "$c: no Allrun"; return 1; }
    t0=$(date +%s)
    (cd "$c" && ./Allrun > log.Allrun 2>&1)
    rc=$?
    echo "$c rc=$rc  $(( $(date +%s) - t0 )) s"
}
export -f run_one

echo "runAll: ${#cases[@]} cases, $jobs at a time"
printf '%s\n' "${cases[@]}" | xargs -P "$jobs" -I{} bash -c 'run_one "$1"' _ {}
python3 getResults.py
