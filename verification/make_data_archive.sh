#!/usr/bin/env bash
# Build the data archive for the non-retracing advection manuscript (Figshare).
#
# Contents:
#   - the repository at HEAD (git archive): library, solver, the three
#     verification workflows, their configs and scripts, the manuscript
#     illustrations (verification/manuscript-illustrations)
#   - every generated result of the three studies (verification/*/results/),
#     which the manuscript's figures and tables cite by path
#   - for every run, the error time series volumeFractionError.dat and the
#     case's system/ directory, so each table entry can be checked without
#     rerunning (solver logs are left out: up to ~10^2 MB per run)
#   - README.md (= verification/DATA_ARCHIVE.md) at the top
#   - MANIFEST.txt: commit, tag, date, host, and a checksum of every file
#
# Run it where the runs are (the Lichtenberg clone), on a tagged commit. It
# refuses to build
#   - from a working tree whose tracked files differ from HEAD,
#   - from an untagged HEAD (ALLOW_UNTAGGED=1 builds a draft for inspection),
#   - while DATA_ARCHIVE.md still carries the DOI placeholder,
#   - if any archived file contains a private path or a credential.
set -euo pipefail
cd "$(dirname "$0")/.."

STUDIES=(reversed-vortex deformation-sphere3d leveque3d)
COMMIT=$(git rev-parse HEAD)
TAG=$(git describe --tags --exact-match HEAD 2>/dev/null || true)
NAME="TwoPhaseFlow-nonretracing-advection-${COMMIT:0:10}"
OUT="${1:-$PWD/$NAME.tar.gz}"

dirty=$(git status --porcelain --untracked-files=no -- . ':(exclude,glob)verification/*/results/**')
if [ -n "$dirty" ]; then
    echo "tracked files differ from HEAD; commit first:" >&2
    echo "$dirty" >&2
    exit 1
fi
if [ -z "$TAG" ]; then
    echo "HEAD carries no tag; tag the release commit first." >&2
    [ "${ALLOW_UNTAGGED:-0}" = 1 ] || exit 1
    TAG="(untagged draft)"
fi
if grep -q '@FIGSHARE_DOI@' verification/DATA_ARCHIVE.md; then
    echo "verification/DATA_ARCHIVE.md still carries the DOI placeholder." >&2
    [ "${ALLOW_UNTAGGED:-0}" = 1 ] || exit 1
fi
for s in "${STUDIES[@]}"; do
    for f in convergence.csv convergence_table.tex courant_check.txt timestep_bound.json; do
        [ -f "verification/$s/results/$f" ] || { echo "missing verification/$s/results/$f" >&2; exit 1; }
    done
done

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
git archive --prefix="$NAME/" HEAD | tar xf - -C "$TMP"
cp verification/DATA_ARCHIVE.md "$TMP/$NAME/README.md"

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
find "$TMP/$NAME" -name __pycache__ -type d -prune -exec rm -rf {} +

# ---- nothing private leaves: paths on the cluster or a workstation, credentials
leaks=$(grep -rIlE "${USER:?}|/work/scratch|/home/[a-z]|/mnt/c/|C:\\\\Users|ghp_[A-Za-z0-9]|github_pat_|PRIVATE KEY" "$TMP/$NAME" || true)
if [ -n "$leaks" ]; then
    echo "private paths or credentials in:" >&2
    echo "$leaks" | sed "s|$TMP/||" >&2
    exit 1
fi

(
    cd "$TMP/$NAME"
    {
        echo "TwoPhaseFlow, non-retracing advection verification: data archive"
        echo "commit   $COMMIT"
        echo "tag      $TAG"
        echo "branch   $(git -C "$OLDPWD" rev-parse --abbrev-ref HEAD)"
        echo "built    $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo
        echo "sha256 of every file:"
        find . -type f ! -name MANIFEST.txt -print0 | sort -z | xargs -0 sha256sum
    } > MANIFEST.txt
)
tar czf "$OUT" -C "$TMP" "$NAME"
(cd "$(dirname "$OUT")" && sha256sum "$(basename "$OUT")" > "$(basename "$OUT").sha256")
echo "wrote $OUT ($(du -h "$OUT" | cut -f1)), commit ${COMMIT:0:10}, tag $TAG"
