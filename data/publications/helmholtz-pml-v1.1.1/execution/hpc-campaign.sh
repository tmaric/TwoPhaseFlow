#!/usr/bin/env bash
source /work/scratch/tm83tomy/ultrasound-cpc-20261115/scripts/hpc-env.sh
set -eo pipefail
campaign="$1"
export CPC_CAMPAIGN="$campaign"
task_out="$CPC_ROOT/results/$campaign"
mkdir -p "$task_out"
test "$(command -v acousticHelmholtzFoam)" = "$FOAM_USER_APPBIN/acousticHelmholtzFoam"
printf 'CAMPAIGN %s JOB %s HOST %s SOURCE %s\n' "$campaign" "$SLURM_JOB_ID" "$(hostname)" "$(git -C "$CPC_REPO" rev-parse HEAD)"
sha256sum "$FOAM_USER_APPBIN/acousticHelmholtzFoam" "$FOAM_USER_LIBBIN/libacousticInterface.so" > "$task_out/binary-sha256.txt"
case "$campaign" in
    homogeneous) template=homogeneousPlaneWave2D ;;
    layered) template=layeredInterface1D ;;
    sphere) template=penetrableSphere ;;
    sphere-pml) template=penetrableSpherePML ;;
    piston) template=pistonRadiation ;;
    gorkov) template=Andrade2019/WedgeLevitatorDropGorkov ;;
    height-*) template=Andrade2019/WedgeLevitatorDropAlpha ;;
    *) echo 'Unknown campaign' >&2; exit 2 ;;
esac
test ! -e "$task_out/case" || { echo 'Refusing to overwrite campaign case' >&2; exit 2; }
mkdir "$task_out/case"
rsync -a "$CPC_REPO/run/acousticTests/FrequencyDomainTests/$template/" "$task_out/case/"
cd "$task_out/case"
case "$campaign" in
    homogeneous) ./run_convergence.sh ;;
    layered)
        python3 run_face_area_comparison.py --output "$task_out/comparison" --ranks 1,4,8
        ;;
    sphere) ./Allrun ;;
    sphere-pml) ./Allrun ;;
    piston) ./meshConv.sh ;;
    gorkov)
        KEEP_CASES=1 ./runMeshConvergence
        KEEP_CASES=1 ./runGorkovStudy
        KEEP_CASES=1 ./runPositionSweep
        ;;
    height-empty) FIG4_EMPTY=1 FIG4_ASPECT_RATIOS='' ./runHeightSweep.sh ;;
    height-1) FIG4_EMPTY=0 FIG4_ASPECT_RATIOS=1 ./runHeightSweep.sh ;;
    height-2) FIG4_EMPTY=0 FIG4_ASPECT_RATIOS=2 ./runHeightSweep.sh ;;
    height-3) FIG4_EMPTY=0 FIG4_ASPECT_RATIOS=3 ./runHeightSweep.sh ;;
esac
printf 'CPC_CAMPAIGN_COMPLETE %s\n' "$campaign" | tee "$task_out/COMPLETED"
