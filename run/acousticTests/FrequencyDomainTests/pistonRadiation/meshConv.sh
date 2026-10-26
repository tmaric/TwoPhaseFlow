#!/bin/sh
set -eu

cd "${0%/*}" || exit 1

resolutions="${MESH_CONV_RESOLUTIONS:-20 30 40 60 80}"
outRoot="${MESH_CONV_OUT:-meshConvergence}"
postArgs="${POSTPROCESS_ARGS:-}"

# The exact reference must pass its checks before any comparison is made.
python3 validate_reference.py

mkdir -p "$outRoot"
printf "cellsPerWavelength,h_over_lambda,time,field_relL2,field_magnitude_relL2,axis_relL2,axis_relLinf,axis_magnitude_relL2\n" \
    > "$outRoot/metrics_summary.csv"

latest_time()
{
    find . -maxdepth 1 -type d -name '[0-9]*' -printf '%f\n' \
        | awk '/^[0-9]+([.][0-9]+)?$/' \
        | sort -g \
        | tail -1
}

metric_value()
{
    key="$1"
    file="$2"
    awk -F= -v key="$key" '
        $1 ~ "^[[:space:]]*" key "[[:space:]]*$" {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2);
            print $2;
            exit
        }
    ' "$file"
}

for n in $resolutions; do
    echo "=== pistonRadiation mesh convergence: ${n} cells per wavelength ==="

    ./Allclean
    rm -rf processor* VTK postProcessing

    MESH_CELLS_PER_LAMBDA="$n" ./Allrun
    if [ -n "$postArgs" ]; then
        # Re-run comparison with caller-specified options.
        # shellcheck disable=SC2086
        python3 postprocess_compare.py $postArgs
    fi

    t="$(latest_time)"
    compareDir="postProcessing/analyticalCompare/$t"
    metricsFile="$compareDir/metrics.txt"

    if [ ! -f "$metricsFile" ]; then
        echo "Error: expected metrics file not found: $metricsFile" >&2
        exit 1
    fi

    caseOut="$outRoot/N${n}"
    rm -rf "$caseOut"
    mkdir -p "$caseOut"
    cp -r "$compareDir" "$caseOut/analyticalCompare"
    cp system/blockMeshDict "$caseOut/"
    cp constant/transportProperties "$caseOut/"
    cp system/decomposeParDict "$caseOut/"
    cp log.* "$caseOut/" 2>/dev/null || true

    hOverLambda="$(awk -v n="$n" 'BEGIN { printf "%.12g", 1.0/n }')"
    row="$n,$hOverLambda,$t"
    for key in field_relL2 field_magnitude_relL2 axis_relL2 axis_relLinf axis_magnitude_relL2; do
        row="$row,$(metric_value "$key" "$metricsFile")"
    done
    echo "$row" >> "$outRoot/metrics_summary.csv"
done

python3 - "$outRoot" $resolutions <<'PY'
import csv
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


out_root = Path(sys.argv[1])
resolutions = sys.argv[2:]

rows = []
for n in resolutions:
    case_dir = out_root / f"N{n}" / "analyticalCompare"
    on_axis = case_dir / "onAxisComparison.csv"
    arc = case_dir / "arcComparison.csv"
    if not on_axis.exists() or not arc.exists():
        raise SystemExit(f"Missing comparison CSVs for N{n}")
    rows.append((n, np.genfromtxt(on_axis, delimiter=",", names=True),
                 np.genfromtxt(arc, delimiter=",", names=True)))

fig, ax = plt.subplots(figsize=(7.4, 4.8))
for i, (n, on, _arc) in enumerate(rows):
    if i == 0:
        ax.plot(on["z_over_rayleigh"], on["p_exact_over_p0"], "k-", lw=2.2, label="Exact")
    ax.plot(on["z_over_rayleigh"], on["p_sim_over_p0"], "--", lw=1.5, label=f"N={n}")
ax.set_xscale("log")
ax.set_xlabel("z / R0")
ax.set_ylabel("|p| / (rho*c*u0)")
ax.grid(True, alpha=0.35)
ax.legend(loc="best")
ax.set_title("pistonRadiation on-axis mesh convergence")
fig.tight_layout()
fig.savefig(out_root / "onAxis_meshConvergence.png", dpi=180)

fig, ax = plt.subplots(figsize=(7.4, 4.8))
for i, (n, _on, arc) in enumerate(rows):
    if i == 0:
        ax.plot(arc["theta_deg"], arc["SPL_exact_dB"], "k-", lw=2.2, label="Exact")
    ax.plot(arc["theta_deg"], arc["SPL_sim_dB"], "--", lw=1.5, label=f"N={n}")
ax.set_xlim(0, 90)
ax.set_xlabel(r"$\theta$ [$^\circ$] from the symmetry axis")
ax.set_ylabel(r"$L_p$ [dB re $20\,\mu$Pa]")
ax.grid(True, alpha=0.35)
ax.legend(loc="best")
ax.set_title(f"Pressure level on the arc r = {rows[0][2]['r_arc_m'][0]:g} m")
fig.tight_layout()
fig.savefig(out_root / "arc_SPL_meshConvergence.png", dpi=180)

with (out_root / "metrics_summary.csv").open(newline="") as f:
    metrics = sorted(csv.DictReader(f), key=lambda row: float(row["cellsPerWavelength"]))

def order(key, i):
    if i == 0:
        return "--"
    e0, e1 = float(metrics[i - 1][key]), float(metrics[i][key])
    n0, n1 = float(metrics[i - 1]["cellsPerWavelength"]), float(metrics[i]["cellsPerWavelength"])
    return f"{math.log(e0 / e1) / math.log(n1 / n0):.2f}"

tex_path = out_root / "piston_convergence_table.tex"
with tex_path.open("w", encoding="utf-8") as f:
    f.write("\\begin{table}[htbp]\n")
    f.write("\\centering\n")
    f.write("\\begin{tabular}{ccccccc}\n")
    f.write("\\hline\n")
    f.write("Cells per wavelength & $h/\\lambda$ & $E_P^{\\Omega}$ & Order"
            " & $E_P^{\\mathrm{axis}}$ & $E_\\infty^{\\mathrm{axis}}$ & Order \\\\\n")
    f.write("\\hline\n")
    for i, row in enumerate(metrics):
        f.write(
            f"{row['cellsPerWavelength']} & "
            f"{float(row['h_over_lambda']):.4f} & "
            f"{float(row['field_relL2']):.3e} & {order('field_relL2', i)} & "
            f"{float(row['axis_relL2']):.3e} & "
            f"{float(row['axis_relLinf']):.3e} & {order('axis_relL2', i)} \\\\\n"
        )
    f.write("\\hline\n")
    f.write("\\end{tabular}\n")
    f.write("\\caption{Mesh convergence of the baffled-piston pressure against the exact"
            " Rayleigh field: complex-pressure errors in the physical region and in the"
            " cell row next to the symmetry axis.}\n")
    f.write("\\label{tab:pistonConvergence}\n")
    f.write("\\end{table}\n")

fig, ax = plt.subplots(figsize=(7.0, 4.6))
x = np.array([float(row["cellsPerWavelength"]) for row in metrics])
for key, mark, label in [
    ("field_relL2", "o-", "physical region"),
    ("axis_relL2", "s-", "axis cell row"),
]:
    ax.loglog(x, [float(row[key]) for row in metrics], mark, label=label)
e0 = float(metrics[0]["field_relL2"])
ax.loglog(x, e0 * (x[0] / x) ** 2, "k:", label="second order")
ax.set_xlabel("cells per wavelength")
ax.set_ylabel("relative complex-pressure error")
ax.grid(True, which="both", alpha=0.35)
ax.legend(loc="best")
ax.set_title("pistonRadiation mesh convergence")
fig.tight_layout()
fig.savefig(out_root / "metrics_meshConvergence.png", dpi=180)
PY

echo "Mesh convergence outputs written to: $outRoot"
echo "Summary: $outRoot/metrics_summary.csv"
