#!/usr/bin/env bash
# Build the data archive for the non-retracing advection manuscript (TUdatalib).
#
# Contents:
#   - the repository at HEAD (git archive): library, solver, the three
#     verification workflows, their configs and scripts
#   - every generated result of the three studies (verification/*/results/),
#     which the manuscript's figures and tables cite by path
#   - for every run, the error time series volumeFractionError.dat and the
#     case's system/ directory, so each table entry can be checked without
#     rerunning (solver logs are left out: up to ~10^2 MB per run)
#   - MANIFEST.txt: commit, date, host, and a checksum of every file
#
# Run it where the runs are (the Lichtenberg clone). It refuses to build from a
# working tree whose tracked files differ from HEAD, because the archive must
# be reproducible from the commit it names.
set -euo pipefail
cd "$(dirname "$0")/.."

STUDIES=(reversed-vortex deformation-sphere3d leveque3d)
COMMIT=$(git rev-parse HEAD)
NAME="TwoPhaseFlow-nonretracing-advection-${COMMIT:0:10}"
OUT="${1:-$PWD/$NAME.tar.gz}"

dirty=$(git status --porcelain --untracked-files=no -- . ':(exclude,glob)verification/*/results/**')
if [ -n "$dirty" ]; then
    echo "tracked files differ from HEAD; commit first:" >&2
    echo "$dirty" >&2
    exit 1
fi
for s in "${STUDIES[@]}"; do
    for f in convergence.csv convergence_table.tex courant_check.txt timestep_bound.json; do
        [ -f "verification/$s/results/$f" ] || { echo "missing verification/$s/results/$f" >&2; exit 1; }
    done
done

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
git archive --prefix="$NAME/" HEAD | tar xf - -C "$TMP"

for s in "${STUDIES[@]}"; do
    mkdir -p "$TMP/$NAME/verification/$s/results"
    cp -r "verification/$s/results/." "$TMP/$NAME/verification/$s/results/"
    for d in verification/$s/runs/*/; do
        c=$(basename "$d")
        [ -f "$d/volumeFractionError.dat" ] || continue
        dst="$TMP/$NAME/verification/$s/runs/$c"
        mkdir -p "$dst"
        cp "$d/volumeFractionError.dat" "$dst/"
        cp -r "$d/system" "$dst/"
    done
done

(
    cd "$TMP/$NAME"
    {
        echo "TwoPhaseFlow, non-retracing advection verification: data archive"
        echo "commit   $COMMIT"
        echo "branch   $(git -C "$OLDPWD" rev-parse --abbrev-ref HEAD)"
        echo "built    $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname)"
        echo
        echo "sha256 of every file:"
        find . -type f ! -name MANIFEST.txt -print0 | sort -z | xargs -0 sha256sum
    } > MANIFEST.txt
)
tar czf "$OUT" -C "$TMP" "$NAME"
echo "wrote $OUT ($(du -h "$OUT" | cut -f1)), commit ${COMMIT:0:10}"
