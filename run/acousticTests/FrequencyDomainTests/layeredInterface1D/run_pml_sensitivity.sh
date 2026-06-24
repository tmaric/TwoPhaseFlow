#!/bin/sh
set -e

case_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
parent_dir="$(dirname -- "$case_dir")"
case_name="$(basename -- "$case_dir")"
out_dir="$case_dir/postProcessing/pmlSensitivity"

sigmas="${PML_SIGMAS:-10000 50000 100000 200000 500000}"

mkdir -p "$out_dir"
rm -f "$out_dir/metrics_summary.csv"
printf 'sigmaMax,Pre_relL2,Pim_relL2,pAbs_relL2\n' > "$out_dir/metrics_summary.csv"

for sigma in $sigmas; do
    tag="$(printf '%s' "$sigma" | sed 's/^-//; s/[^0-9A-Za-z_]/_/g')"
    work_dir="$parent_dir/${case_name}_sigma_${tag}"

    rm -rf "$work_dir"
    mkdir -p "$work_dir"
    cp -a "$case_dir/0.orig" "$case_dir/constant" "$case_dir/system" "$work_dir/"
    cp -a "$case_dir/Allrun" "$case_dir/Allclean" "$case_dir/prepareCase" \
        "$case_dir/caseParams.sh" "$case_dir/postprocess_interface.py" "$work_dir/"

    python3 - "$work_dir/caseParams.sh" "$sigma" <<'PY'
import pathlib
import re
import sys

path = pathlib.Path(sys.argv[1])
sigma = sys.argv[2]
text = path.read_text(encoding="utf-8")
text = re.sub(r"(?m)^SIGMA_MAX=.*$", f"SIGMA_MAX={sigma}", text)
path.write_text(text, encoding="utf-8")
PY

    echo "Running PML sensitivity case sigmaMax=$sigma"
    (cd "$work_dir" && ./Allrun)

    metrics="$work_dir/postProcessing/interfaceValidation/1/metrics.txt"
    pre="$(awk -F ' = ' '$1=="Pre_relL2"{print $2}' "$metrics")"
    pim="$(awk -F ' = ' '$1=="Pim_relL2"{print $2}' "$metrics")"
    pabs="$(awk -F ' = ' '$1=="|p|_relL2"{print $2}' "$metrics")"
    printf '%s,%s,%s,%s\n' "$sigma" "$pre" "$pim" "$pabs" >> "$out_dir/metrics_summary.csv"

    mkdir -p "$out_dir/sigma_${tag}"
    cp -a "$work_dir/postProcessing/interfaceValidation/1/." "$out_dir/sigma_${tag}/"

    if [ "${KEEP_PML_SWEEP_CASES:-0}" != "1" ]; then
        rm -rf "$work_dir"
    fi
done

python3 "$case_dir/postprocess_pml_sensitivity.py"

echo "Wrote $out_dir/metrics_summary.csv"
echo "Wrote $out_dir/pml_sigma_sensitivity.png"
