#!/bin/sh
set -eu
cd "${0%/*}" || exit 1

resolutions="${BASELINE_RESOLUTIONS:-16 24 32 48 64 96}"
families="${BASELINE_MESH_FAMILIES:-orthogonal warpedInterior}"
boundaries="${BASELINE_BOUNDARY_MODES:-dirichlet mixed}"
out_dir=studyResults/convergence
mkdir -p "$out_dir"
printf '%s\n' 'meshFamily,boundaryMode,cellsPerWavelength,cells,hOverLambda,maxNonOrthogonality,avgNonOrthogonality,pressureRelL2,pressureRelLinf,velocityRelL2,velocityInteriorRelL2,velocityBoundaryRelL2' > "$out_dir/metrics_summary.csv"

for boundary in $boundaries; do
    for family in $families; do
        for resolution in $resolutions; do
            echo "=== baseline: $boundary, $family, Nlambda=$resolution ==="
            CELLS_PER_WAVELENGTH=$resolution MESH_FAMILY=$family BOUNDARY_MODE=$boundary ./Allrun
            metrics=postProcessing/homogeneousBaseline/1/metrics.csv
            tail -n 1 "$metrics" | awk -v family="$family" -v boundary="$boundary" 'BEGIN{FS=OFS=","} {print family,boundary,$0}' >> "$out_dir/metrics_summary.csv"
        done
    done
done

python3 summarize_convergence.py "$out_dir/metrics_summary.csv" "$out_dir"
