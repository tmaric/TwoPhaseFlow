#!/bin/sh
set -eu

cd "${0%/*}" || exit 1
. ./caseParams.sh

resultsDir="sweepResults"
mkdir -p "${resultsDir}"

csvFile="${resultsDir}/heightSweepResults.csv"
printf "index,height_fac,D_m,Fn_reflector1_N\n" > "${csvFile}"

lambda="$(awk -v c="$SOUND_SPEED" -v f="$DRIVE_F" 'BEGIN{printf "%.16g", c/f}')"

if [ "$SWEEP_POINTS" -lt 2 ]; then
    echo "SWEEP_POINTS must be >= 2" >&2
    exit 1
fi

echo "Running height sweep: ${SWEEP_POINTS} points from ${SWEEP_FAC_MIN} to ${SWEEP_FAC_MAX}"

i=0
while [ "$i" -lt "$SWEEP_POINTS" ]; do
    fac="$(awk -v i="$i" -v n="$SWEEP_POINTS" -v a="$SWEEP_FAC_MIN" -v b="$SWEEP_FAC_MAX" 'BEGIN{printf "%.16g", a + (b-a)*i/(n-1)}')"
    dVal="$(awk -v ff="$fac" -v lam="$lambda" 'BEGIN{printf "%.16g", ff*lam}')"

    echo "[$((i+1))/${SWEEP_POINTS}] HEIGHT_FAC=${fac} D=${dVal} m"
    HEIGHT_FAC="$fac" ./Allrun > "${resultsDir}/log_height_$(printf "%03d" "$i").txt" 2>&1

    forceFile="postProcessing/reflectorForce/force.dat"
    if [ ! -f "${forceFile}" ]; then
        echo "Missing ${forceFile} for HEIGHT_FAC=${fac}" >&2
        exit 1
    fi

    fn="$(awk 'BEGIN{v=""} /^[[:space:]]*#/ {next} NF>=2 {v=$2} END{if(v==""){exit 1}else{print v}}' "${forceFile}")"
    printf "%d,%.16g,%.16g,%.16g\n" "$i" "$fac" "$dVal" "$fn" >> "${csvFile}"

    i=$((i+1))
done

python3 plotHeightSweep.py "${csvFile}" "${resultsDir}/heightSweep_force.png"
echo "Sweep complete. Results: ${csvFile}"
echo "Plot written to: ${resultsDir}/heightSweep_force.png"
