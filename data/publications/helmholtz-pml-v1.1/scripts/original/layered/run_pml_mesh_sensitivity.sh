#!/bin/sh
set -e

case_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
parent_dir="$(dirname -- "$case_dir")"
case_name="$(basename -- "$case_dir")"
out_dir="$case_dir/postProcessing/pmlMeshSensitivity"
. "$case_dir/caseParams.sh"

sigma="${PML_MESH_SIGMA:-500000}"
nx_values="${PML_MESH_NX_VALUES:-560 1120 2240 4480 8960}"

mkdir -p "$out_dir"
rm -f "$out_dir/metrics_summary.csv"
printf 'sigmaMax,NX,dx_m,P_relL2,Pre_relL2,Pim_relL2\n' > "$out_dir/metrics_summary.csv"

for nx in $nx_values; do
    tag="$(printf '%s' "$nx" | sed 's/[^0-9A-Za-z_]/_/g')"
    work_dir="$parent_dir/${case_name}_sigmaMesh_NX_${tag}"

    rm -rf "$work_dir"
    mkdir -p "$work_dir"
    cp -a "$case_dir/0.orig" "$case_dir/0.templates" \
        "$case_dir/constant" "$case_dir/system" "$work_dir/"
    cp -a "$case_dir/Allrun" "$case_dir/Allclean" "$case_dir/prepareCase" \
        "$case_dir/caseParams.sh" "$case_dir/postprocess_interface.py" "$work_dir/"

    python3 - "$work_dir/caseParams.sh" "$sigma" "$nx" <<'PY'
import pathlib
import re
import sys

path = pathlib.Path(sys.argv[1])
sigma = sys.argv[2]
nx = sys.argv[3]
text = path.read_text(encoding="utf-8")
text = re.sub(r"(?m)^SIGMA_MAX=.*$", f"SIGMA_MAX={sigma}", text)
text = re.sub(r"(?m)^NX=.*$", f"NX={nx}", text)
path.write_text(text, encoding="utf-8")
PY

    echo "Running PML mesh sensitivity case sigmaMax=$sigma NX=$nx"
    (cd "$work_dir" && ./Allrun)

    metrics="$work_dir/postProcessing/interfaceValidation/1/metrics.txt"
    pressure="$(awk -F ' = ' '$1=="P_relL2"{print $2}' "$metrics")"
    pre="$(awk -F ' = ' '$1=="Pre_relL2"{print $2}' "$metrics")"
    pim="$(awk -F ' = ' '$1=="Pim_relL2"{print $2}' "$metrics")"
    dx="$(awk -v xmax="$X_MAX" -v xmin="$X_MIN" -v n="$nx" \
        'BEGIN{printf "%.12e", (xmax-xmin)/n}')"
    printf '%s,%s,%s,%s,%s,%s\n' \
        "$sigma" "$nx" "$dx" "$pressure" "$pre" "$pim" \
        >> "$out_dir/metrics_summary.csv"

    mkdir -p "$out_dir/NX_${tag}"
    cp -a "$work_dir/postProcessing/interfaceValidation/1/." "$out_dir/NX_${tag}/"

    if [ "${KEEP_PML_MESH_SWEEP_CASES:-0}" != "1" ]; then
        rm -rf "$work_dir"
    fi
done

python3 "$case_dir/postprocess_pml_mesh_sensitivity.py"

echo "Wrote $out_dir/metrics_summary.csv"
echo "Wrote $out_dir/pml_mesh_sensitivity.png"
