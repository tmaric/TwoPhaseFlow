#!/bin/sh
set -e

case_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
parent_dir="$(dirname -- "$case_dir")"
case_name="$(basename -- "$case_dir")"
out_dir="$case_dir/studyResults/interfaceConvergence"
nx_values="${INTERFACE_NX_VALUES:-160 240 320 480 640 960}"
offsets="${INTERFACE_OFFSETS:-0 0.25 0.5}"
xmin=0
xmax=0.14

mkdir -p "$out_dir"
printf '%s\n' 'offset,NX,dx,interfaceX,pressureRelL2,RReal,RImag,RAbsError,TReal,TImag,TAbsError' > "$out_dir/metrics_summary.csv"

for offset in $offsets; do
    for nx in $nx_values; do
        dx=$(awk -v a="$xmin" -v b="$xmax" -v n="$nx" 'BEGIN{printf "%.16g",(b-a)/n}')
        xi=$(awk -v a="$xmin" -v b="$xmax" -v d="$dx" -v o="$offset" 'BEGIN{printf "%.16g",0.5*(a+b)+o*d}')
        tag="o${offset}_N${nx}"
        work_dir="$parent_dir/${case_name}_interface_${tag}"
        rm -rf "$work_dir"
        mkdir -p "$work_dir"
        cp -a "$case_dir/0.orig" "$case_dir/0.templates" "$case_dir/constant" "$case_dir/system" "$work_dir/"
        cp -a "$case_dir/Allrun" "$case_dir/Allclean" "$case_dir/prepareCase" \
            "$case_dir/caseParams.sh" "$case_dir/postprocess_interface.py" \
            "$case_dir/postprocess_interface_convergence.py" "$work_dir/"
        echo "Running interface convergence offset=$offset NX=$nx"
        (cd "$work_dir" && CASE_MODE=interface X_MIN=$xmin X_MAX=$xmax PML_XMAX=$xmax \
            X_INTERFACE=$xi NX=$nx SIGMA_MAX=0 ./Allrun)
        metrics="$work_dir/postProcessing/interfaceConvergence/1/metrics.csv"
        tail -n 1 "$metrics" | awk -v offset="$offset" 'BEGIN{FS=OFS=","}{print offset,$0}' >> "$out_dir/metrics_summary.csv"
        rm -rf "$work_dir"
    done
done

python3 "$case_dir/summarize_interface_convergence.py" "$out_dir/metrics_summary.csv" "$out_dir"
