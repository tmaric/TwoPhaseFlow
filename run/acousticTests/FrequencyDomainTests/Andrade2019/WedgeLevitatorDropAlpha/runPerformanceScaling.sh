#!/usr/bin/env bash
set -eo pipefail

cd "$(dirname "$0")"
. "${WM_PROJECT_DIR:?}/bin/tools/RunFunctions"
. ./caseParams.sh

# Performance studies:
#   strong: fixed mesh, varying MPI ranks
#   size:   fixed MPI ranks, varying cfMesh maximum cell size
#   all:    run both studies
#
# Examples:
#   ./runPerformanceScaling.sh
#   STRONG_RANKS="1 2 4 8 16" ./runPerformanceScaling.sh strong
#   SIZE_CELL_SIZES="8e-5 6e-5 4e-5 3e-5" SIZE_RANKS=8 ./runPerformanceScaling.sh size

MODE="${1:-all}"
STRONG_RANKS="${STRONG_RANKS:-1 2 4 8}"
STRONG_CELL_SIZE="${STRONG_CELL_SIZE:-${CFMESH_MAX_CELL_SIZE}}"
SIZE_CELL_SIZES="${SIZE_CELL_SIZES:-8e-5 6e-5 4e-5 3e-5}"
SIZE_RANKS="${SIZE_RANKS:-4}"
MPIEXEC="${MPIEXEC:-mpirun}"
MPI_OPTIONS="${MPI_OPTIONS:-}"
RESULTS_DIR="${PERFORMANCE_RESULTS_DIR:-performanceResults/$(date +%Y%m%d-%H%M%S)}"
PERFORMANCE_WARMUP="${PERFORMANCE_WARMUP:-1}"
WARMED_UP=0

if [ -z "$MPI_OPTIONS" ] && [[ "${FOAM_MPI:-}" == *openmpi* ]]; then
    MPI_OPTIONS="--oversubscribe"
fi
read -r -a MPI_OPTION_ARRAY <<< "$MPI_OPTIONS"

case "$MODE" in
    strong|size|all) ;;
    *)
        echo "Usage: $0 [strong|size|all]" >&2
        exit 2
        ;;
esac

need_cmd()
{
    command -v "$1" >/dev/null 2>&1 || {
        echo "Error: required command '$1' not found in PATH" >&2
        exit 1
    }
}

for commandName in \
    cartesian2DMesh extrudeMesh setAlphaField setPMLFields \
    decomposePar checkMesh foamDictionary acousticHelmholtzFoam \
    "$MPIEXEC" /usr/bin/time
do
    need_cmd "$commandName"
done

mkdir -p "$RESULTS_DIR/logs"
RESULTS_DIR="$(cd "$RESULTS_DIR" && pwd)"

decomposeBackup="$(mktemp)"
cp system/decomposeParDict "$decomposeBackup"
restore_decompose_dict()
{
    cp "$decomposeBackup" system/decomposeParDict
    rm -f "$decomposeBackup"
}
trap restore_decompose_dict EXIT

prepare_mesh()
{
    local cellSize="$1"
    local tag="$2"

    echo "Preparing mesh '$tag' with CFMESH_MAX_CELL_SIZE=$cellSize"
    ./Allclean
    CFMESH_MAX_CELL_SIZE="$cellSize" ./prepareCase

    cartesian2DMesh > "$RESULTS_DIR/logs/${tag}.cartesian2DMesh.log" 2>&1
    extrudeMesh > "$RESULTS_DIR/logs/${tag}.extrudeMesh.log" 2>&1
    ./postMesh

    restore0Dir
    setAlphaField > "$RESULTS_DIR/logs/${tag}.setAlphaField.log" 2>&1
    CFMESH_MAX_CELL_SIZE="$cellSize" ./refineAlphaInterface \
        > "$RESULTS_DIR/logs/${tag}.refineAlphaInterface.log" 2>&1
    rm -f 0/Pim.in
    setPMLFields > "$RESULTS_DIR/logs/${tag}.setPMLFields.log" 2>&1

    checkMesh -constant > "$RESULTS_DIR/logs/${tag}.checkMesh.log" 2>&1
}

mesh_cells()
{
    awk '/^[[:space:]]*cells:/ {print $2; exit}' "$1"
}

solver_clock()
{
    awk '
        /ExecutionTime =/ {
            for (i = 1; i <= NF; ++i)
            {
                if ($i == "ClockTime")
                {
                    clock = $(i + 2)
                }
            }
        }
        END {print clock}
    ' "$1"
}

run_solver()
{
    local ranks="$1"
    local tag="$2"
    local logFile="$RESULTS_DIR/logs/${tag}.solver.log"
    local timeFile="$RESULTS_DIR/logs/${tag}.elapsed"

    rm -rf processor* 1 postProcessing
    foamDictionary system/decomposeParDict \
        -entry numberOfSubdomains -set "$ranks" >/dev/null

    if [ "$ranks" -eq 1 ]; then
        /usr/bin/time -f '%e' -o "$timeFile" \
            acousticHelmholtzFoam > "$logFile" 2>&1
    else
        decomposePar -force > "$RESULTS_DIR/logs/${tag}.decomposePar.log" 2>&1
        /usr/bin/time -f '%e' -o "$timeFile" \
            "$MPIEXEC" "${MPI_OPTION_ARRAY[@]}" -np "$ranks" \
            acousticHelmholtzFoam -parallel \
            > "$logFile" 2>&1
    fi

    if ! grep -q '^End$' "$logFile"; then
        echo "Error: solver did not finish normally; see $logFile" >&2
        exit 1
    fi

    RUN_ELAPSED="$(cat "$timeFile")"
    RUN_SOLVER_CLOCK="$(solver_clock "$logFile")"
}

warm_up_once()
{
    local ranks="$1"
    local tag="$2"

    if [ "$PERFORMANCE_WARMUP" = "0" ] || [ "$WARMED_UP" = "1" ]; then
        return
    fi

    echo "Warm-up run: ranks=$ranks"
    run_solver "$ranks" "${tag}.warmup"
    WARMED_UP=1
}

run_strong_scaling()
{
    local csv="$RESULTS_DIR/strong_scaling.csv"
    local checkLog="$RESULTS_DIR/logs/strong.checkMesh.log"
    local cells unknowns ranks tag speedup efficiency
    local baselineElapsed=""
    local baselineRanks=""

    prepare_mesh "$STRONG_CELL_SIZE" strong
    cells="$(mesh_cells "$checkLog")"
    unknowns=$((2*cells))
    set -- $STRONG_RANKS
    warm_up_once "$1" strong

    printf '%s\n' \
        'cell_size,cells,unknowns,ranks,elapsed_s,solver_clock_s,speedup,efficiency' \
        > "$csv"

    for ranks in $STRONG_RANKS; do
        tag="strong.np${ranks}"
        echo "Strong scaling: cells=$cells ranks=$ranks"
        run_solver "$ranks" "$tag"

        if [ -z "$baselineElapsed" ]; then
            baselineElapsed="$RUN_ELAPSED"
            baselineRanks="$ranks"
        fi
        speedup="$(awk -v base="$baselineElapsed" -v elapsed="$RUN_ELAPSED" \
            'BEGIN {printf "%.6g", base/elapsed}')"
        efficiency="$(awk -v speedup="$speedup" -v baseRanks="$baselineRanks" -v ranks="$ranks" \
            'BEGIN {printf "%.6g", speedup*baseRanks/ranks}')"

        printf '%s,%s,%s,%s,%s,%s,%s,%s\n' \
            "$STRONG_CELL_SIZE" "$cells" "$unknowns" "$ranks" \
            "$RUN_ELAPSED" "$RUN_SOLVER_CLOCK" "$speedup" "$efficiency" \
            >> "$csv"
    done
}

run_size_scaling()
{
    local csv="$RESULTS_DIR/problem_size_scaling.csv"
    local cellSize safeSize tag checkLog cells unknowns

    printf '%s\n' \
        'cell_size,cells,unknowns,ranks,elapsed_s,solver_clock_s' \
        > "$csv"

    for cellSize in $SIZE_CELL_SIZES; do
        safeSize="$(printf '%s' "$cellSize" | tr '.+-' 'p__')"
        tag="size.${safeSize}"
        prepare_mesh "$cellSize" "$tag"

        checkLog="$RESULTS_DIR/logs/${tag}.checkMesh.log"
        cells="$(mesh_cells "$checkLog")"
        unknowns=$((2*cells))

        echo "Problem-size scaling: cellSize=$cellSize cells=$cells ranks=$SIZE_RANKS"
        warm_up_once "$SIZE_RANKS" "$tag"
        run_solver "$SIZE_RANKS" "${tag}.np${SIZE_RANKS}"
        printf '%s,%s,%s,%s,%s,%s\n' \
            "$cellSize" "$cells" "$unknowns" "$SIZE_RANKS" \
            "$RUN_ELAPSED" "$RUN_SOLVER_CLOCK" >> "$csv"
    done
}

case "$MODE" in
    strong) run_strong_scaling ;;
    size) run_size_scaling ;;
    all)
        run_strong_scaling
        run_size_scaling
        ;;
esac

echo "Performance results written to $RESULTS_DIR"
