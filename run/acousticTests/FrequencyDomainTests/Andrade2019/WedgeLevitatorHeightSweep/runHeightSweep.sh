#!/bin/sh
set -eu

cd "${0%/*}" || exit 1
. ./caseParams.sh

resultsDir="${RESULTS_DIR}"
mkdir -p "${resultsDir}"

csvFile="${resultsDir}/heightSweepResults.csv"
fig4Base="${resultsDir}/${FIG4_OUTPUT_BASENAME}"
printf "series_key,series_label,drop_mode,aspect_ratio,index,height_mm,D_m,height_fac,Fn_reflector1_wedge_N,Fn_reflector1_axisym_N,axisym_factor\n" > "${csvFile}"

lambda="$(awk -v c="$SOUND_SPEED" -v f="$DRIVE_F" 'BEGIN{printf "%.16g", c/f}')"
axisymFactor="$(awk -v deg="$WEDGE_DEG" 'BEGIN{if (deg <= 0) {exit 1} printf "%.16g", 360.0/deg}')"

if [ "$SWEEP_POINTS" -lt 2 ]; then
    echo "SWEEP_POINTS must be >= 2" >&2
    exit 1
fi

echo "Running height sweep: ${SWEEP_POINTS} points from ${SWEEP_FAC_MIN} to ${SWEEP_FAC_MAX}"

run_series()
{
    seriesKey="$1"
    seriesLabel="$2"
    dropMode="$3"
    aspectRatio="$4"

    echo "Series ${seriesLabel}"

    i=0
    while [ "$i" -lt "$SWEEP_POINTS" ]; do
        fac="$(awk -v i="$i" -v n="$SWEEP_POINTS" -v a="$SWEEP_FAC_MIN" -v b="$SWEEP_FAC_MAX" 'BEGIN{printf "%.16g", a + (b-a)*i/(n-1)}')"
        dVal="$(awk -v ff="$fac" -v lam="$lambda" 'BEGIN{printf "%.16g", ff*lam}')"
        dMm="$(awk -v d="$dVal" 'BEGIN{printf "%.16g", 1000.0*d}')"

        echo "[$((i+1))/${SWEEP_POINTS}] ${seriesKey} HEIGHT_FAC=${fac} D=${dVal} m (${dMm} mm)"

        # Force a fresh run for each sweep point; RunFunctions otherwise skips
        # applications when log.<app> exists from previous point.
        rm -f \
            log.gmshToFoam \
            log.changeDictionary \
            log.setPMLFields \
            log.setExprFields \
            log.decomposePar \
            log.reconstructPar \
            "log.${SOLVER}" \
            log.acousticHelmholtzFoam
        rm -rf processor*
        rm -rf postProcessing/reflectorForce

        HEIGHT_FAC="$fac" \
        DROP_MODE="$dropMode" \
        DROP_ASPECT_RATIO="$aspectRatio" \
        TPF_INTERNAL_SWEEP_RUN=1 \
        ./Allrun > "${resultsDir}/log_${seriesKey}_height_$(printf "%03d" "$i").txt" 2>&1

        forceFile="postProcessing/reflectorForce/force.dat"
        if [ ! -f "${forceFile}" ]; then
            echo "Missing ${forceFile} for series=${seriesKey} HEIGHT_FAC=${fac}" >&2
            exit 1
        fi

        fnWedge="$(awk 'BEGIN{v=""} /^[[:space:]]*#/ {next} NF>=2 {v=$2} END{if(v==""){exit 1}else{print v}}' "${forceFile}")"
        fnAxisym="$(awk -v f="$fnWedge" -v m="$axisymFactor" 'BEGIN{printf "%.16g", f*m}')"
        printf "%s,%s,%s,%.16g,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n" \
            "$seriesKey" "$seriesLabel" "$dropMode" "$aspectRatio" "$i" \
            "$dMm" "$dVal" "$fac" "$fnWedge" "$fnAxisym" "$axisymFactor" >> "${csvFile}"

        i=$((i+1))
    done
}

run_series "empty" "Empty levitator" "empty" "0"

for aspectRatio in $FIG4_ASPECT_RATIOS; do
    if ! awk -v ar="$aspectRatio" 'BEGIN{exit !(ar > 0)}'; then
        echo "Invalid aspect ratio '${aspectRatio}' in FIG4_ASPECT_RATIOS" >&2
        exit 1
    fi

    seriesKey="ab$(printf "%s" "$aspectRatio" | tr '. ' '__')"
    if awk -v ar="$aspectRatio" 'BEGIN{exit !(ar == 1)}'; then
        seriesLabel="Sphere (a/b=1)"
    else
        seriesLabel="Oblate spheroid (a/b=${aspectRatio})"
    fi

    run_series "$seriesKey" "$seriesLabel" "oblate" "$aspectRatio"
done

python3 plotHeightSweep.py "${csvFile}" "${fig4Base}"
echo "Sweep complete. Results: ${csvFile}"
echo "Fig. 4 reproduction written to: ${fig4Base}.png"
echo "Fig. 4 reproduction written to: ${fig4Base}.pdf"
