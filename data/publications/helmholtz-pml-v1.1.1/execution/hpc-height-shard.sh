#!/usr/bin/env bash
source /work/scratch/tm83tomy/ultrasound-cpc-20261115/scripts/hpc-env.sh
set -eo pipefail
series="$1"; part="$2"; count="$3"; fac_min="$4"; fac_max="$5"
task_out="$CPC_ROOT/results/height-shards/$series-$part"
test ! -e "$task_out" || { echo 'Refusing to overwrite shard' >&2; exit 2; }
mkdir -p "$task_out/case"
rsync -a "$CPC_REPO/run/acousticTests/FrequencyDomainTests/Andrade2019/WedgeLevitatorDropAlpha/" "$task_out/case/"
sha256sum "$FOAM_USER_APPBIN/acousticHelmholtzFoam" "$FOAM_USER_LIBBIN/libacousticInterface.so" > "$task_out/binary-sha256.txt"
printf 'SHARD %s %s JOB %s HOST %s SOURCE %s\n' "$series" "$part" "$SLURM_JOB_ID" "$(hostname)" "$(git -C "$CPC_REPO" rev-parse HEAD)"
export SWEEP_POINTS="$count" SWEEP_FAC_MIN="$fac_min" SWEEP_FAC_MAX="$fac_max"
if [ "$series" = empty ]; then
 export FIG4_EMPTY=1 FIG4_ASPECT_RATIOS=''
else
 export FIG4_EMPTY=0 FIG4_ASPECT_RATIOS="$series"
fi
cd "$task_out/case"
./runHeightSweep.sh
printf 'CPC_HEIGHT_SHARD_COMPLETE %s %s\n' "$series" "$part" | tee "$task_out/COMPLETED"
