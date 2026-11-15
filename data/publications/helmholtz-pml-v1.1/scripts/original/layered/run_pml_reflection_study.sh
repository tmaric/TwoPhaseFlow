#!/bin/sh
set -e

case_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
parent_dir="$(dirname -- "$case_dir")"
case_name="$(basename -- "$case_dir")"
out_dir="$case_dir/studyResults/pmlReflection"
frequency=20000
c=343
rho=1.2
lambda=$(awk -v c="$c" -v f="$frequency" 'BEGIN{printf "%.16g",c/f}')
physical_length=$(awk -v l="$lambda" 'BEGIN{printf "%.16g",5*l}')

mkdir -p "$out_dir"
printf '%s\n' 'study,value,NX,cellsPerWavelength,sigmaOverOmega,pmlThicknessOverLambda,pmlOrder,reflectionMagnitude,fitResidual,outerToInnerAmplitude' > "$out_dir/metrics_summary.csv"

run_case()
{
    study=$1
    value=$2
    nlambda=$3
    sigma_ratio=$4
    thickness_ratio=$5
    order=$6
    pml_length=$(awk -v r="$thickness_ratio" -v l="$lambda" 'BEGIN{printf "%.16g",r*l}')
    xmax=$(awk -v a="$physical_length" -v b="$pml_length" 'BEGIN{printf "%.16g",a+b}')
    nx=$(awk -v xmax="$xmax" -v l="$lambda" -v n="$nlambda" 'BEGIN{printf "%d",int(xmax/l*n+0.5)}')
    sigma=$(awk -v r="$sigma_ratio" -v f="$frequency" 'BEGIN{printf "%.16g",r*2*atan2(0,-1)*f}')
    work_dir="$parent_dir/${case_name}_pml_${study}_${value}"
    rm -rf "$work_dir"
    mkdir -p "$work_dir"
    cp -a "$case_dir/0.orig" "$case_dir/0.templates" "$case_dir/constant" "$case_dir/system" "$work_dir/"
    cp -a "$case_dir/Allrun" "$case_dir/Allclean" "$case_dir/prepareCase" \
        "$case_dir/caseParams.sh" "$case_dir/postprocess_interface.py" \
        "$case_dir/postprocess_interface_convergence.py" "$case_dir/postprocess_pml_reflection.py" "$work_dir/"
    echo "Running PML reflection $study=$value"
    (cd "$work_dir" && CASE_MODE=pmlReflection DRIVE_F=$frequency RHOG=$rho RHOL=$rho CG=$c CL=$c \
        X_MIN=0 X_MAX=$xmax PML_XMAX=$xmax PML_L=$pml_length X_INTERFACE=$physical_length \
        NX=$nx SIGMA_MAX=$sigma PO=$order ./Allrun)
    metrics="$work_dir/postProcessing/pmlReflection/1/metrics.csv"
    tail -n 1 "$metrics" | awk -v study="$study" -v value="$value" 'BEGIN{FS=OFS=","}{print study,value,$0}' >> "$out_dir/metrics_summary.csv"
    rm -rf "$work_dir"
}

for value in 0.5 1 2 4 8; do run_case damping "$value" 40 "$value" 1 3; done
for value in 0.5 0.75 1 1.5; do run_case thickness "$value" 40 4 "$value" 3; done
for value in 20 30 40 60 80; do run_case resolution "$value" "$value" 4 1 3; done
for value in 2 3 4; do run_case order "$value" 40 4 1 "$value"; done

python3 "$case_dir/summarize_pml_reflection.py" "$out_dir/metrics_summary.csv" "$out_dir"
