#!/bin/sh
set -eu

cd "${0%/*}" || exit 1

resolutions="${MESH_CONV_RESOLUTIONS:-20 50 100}"
outRoot="${MESH_CONV_OUT:-meshConvergence}"
postArgs="${POSTPROCESS_ARGS:-}"

mkdir -p "$outRoot"
: > "$outRoot/metrics_summary.csv"
printf "cellsPerWavelength,h_over_lambda,time,relL2,relLinf,absLinf,farField_relL2,farField_relLinf,farField_absLinf\n" \
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
        # Re-run comparison with caller-specified sampling/far-field options.
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
    cp constant/pistonRadiation.geo "$caseOut/"
    cp constant/transportProperties "$caseOut/"
    cp system/decomposeParDict "$caseOut/"
    cp log.* "$caseOut/" 2>/dev/null || true

    relL2="$(metric_value relL2 "$metricsFile")"
    relLinf="$(metric_value relLinf "$metricsFile")"
    absLinf="$(metric_value absLinf "$metricsFile")"
    ffRelL2="$(metric_value farField_relL2 "$metricsFile")"
    ffRelLinf="$(metric_value farField_relLinf "$metricsFile")"
    ffAbsLinf="$(metric_value farField_absLinf "$metricsFile")"
    hOverLambda="$(awk -v n="$n" 'BEGIN { printf "%.12g", 1.0/n }')"
    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
        "$n" "$hOverLambda" "$t" "$relL2" "$relLinf" "$absLinf" "$ffRelL2" "$ffRelLinf" "$ffAbsLinf" \
        >> "$outRoot/metrics_summary.csv"
done

python3 - "$outRoot" $resolutions <<'PY'
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


out_root = Path(sys.argv[1])
resolutions = sys.argv[2:]

summary_rows = []
for n in resolutions:
    case_dir = out_root / f"N{n}" / "analyticalCompare"
    on_axis = case_dir / "onAxisComparison.csv"
    far_field = case_dir / "farFieldPatternComparison.csv"
    if not on_axis.exists() or not far_field.exists():
        raise SystemExit(f"Missing comparison CSVs for N{n}")

    on = np.genfromtxt(on_axis, delimiter=",", names=True)
    ff = np.genfromtxt(far_field, delimiter=",", names=True)
    summary_rows.append((n, on, ff))

fig, ax = plt.subplots(figsize=(7.4, 4.8))
first_on = summary_rows[0][1]
ax.plot(
    np.clip(first_on["z_over_rayleigh"], 1e-12, None),
    first_on["p_analytic_over_p0"],
    "k-",
    lw=2.2,
    label="Analytical",
)
for n, on, _ff in summary_rows:
    ax.plot(
        np.clip(on["z_over_rayleigh"], 1e-12, None),
        on["p_sim_over_p0"],
        "--",
        lw=1.5,
        label=f"N={n}",
    )
ax.set_xscale("log")
ax.set_xlim(5e-4, 1)
ax.set_xlabel("z / R0")
ax.set_ylabel("|p| / (rho*c*u0)")
ax.grid(True, alpha=0.35)
ax.legend(loc="best")
ax.set_title("pistonRadiation near-field/on-axis mesh convergence")
fig.tight_layout()
fig.savefig(out_root / "nearField_onAxis_meshConvergence.png", dpi=180)

fig, ax = plt.subplots(figsize=(7.4, 4.8))
first_ff = summary_rows[0][2]
ax.plot(
    first_ff["theta_deg"],
    first_ff["SPL_analytic_dB"],
    "k-",
    lw=2.2,
    label="Analytical",
)
for n, _on, ff in summary_rows:
    ax.plot(
        ff["theta_deg"],
        ff["SPL_sim_dB"],
        "--",
        lw=1.5,
        label=f"N={n}",
    )
ax.set_xlim(0, 90)
ax.set_xlabel("theta [deg] (from axis)")
ax.set_ylabel("SPL [dB re 20uPa]")
ax.grid(True, alpha=0.35)
ax.legend(loc="best")
ax.set_title("pistonRadiation far-field SPL mesh convergence")
fig.tight_layout()
fig.savefig(out_root / "farField_SPL_meshConvergence.png", dpi=180)

fig, ax = plt.subplots(figsize=(7.4, 4.8))
ax.plot(
    first_ff["theta_deg"],
    first_ff["directivity_analytic"],
    "k-",
    lw=2.2,
    label="Analytical",
)
for n, _on, ff in summary_rows:
    ax.plot(
        ff["theta_deg"],
        ff["directivity_sim"],
        "--",
        lw=1.5,
        label=f"N={n}",
    )
ax.set_xlim(0, 90)
ax.set_xlabel("theta [deg] (from axis)")
ax.set_ylabel("normalized directivity")
ax.grid(True, alpha=0.35)
ax.legend(loc="best")
ax.set_title("pistonRadiation far-field directivity mesh convergence")
fig.tight_layout()
fig.savefig(out_root / "farField_directivity_meshConvergence.png", dpi=180)

metrics_path = out_root / "metrics_summary.csv"
with metrics_path.open(newline="") as f:
    metrics = list(csv.DictReader(f))

metrics.sort(key=lambda row: float(row["cellsPerWavelength"]))

def value(row, key):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")

def observed_orders(rows, key):
    orders = [float("nan")]
    for prev, cur in zip(rows[:-1], rows[1:]):
        e_prev = value(prev, key)
        e_cur = value(cur, key)
        n_prev = value(prev, "cellsPerWavelength")
        n_cur = value(cur, "cellsPerWavelength")
        if e_prev > 0.0 and e_cur > 0.0 and n_prev > 0.0 and n_cur > n_prev:
            orders.append(float(np.log(e_prev / e_cur) / np.log(n_cur / n_prev)))
        else:
            orders.append(float("nan"))
    return orders

order_rel_l2 = observed_orders(metrics, "relL2")
order_rel_linf = observed_orders(metrics, "relLinf")
order_ff_rel_l2 = observed_orders(metrics, "farField_relL2")
order_ff_rel_linf = observed_orders(metrics, "farField_relLinf")

for i, row in enumerate(metrics):
    row["order_relL2"] = order_rel_l2[i]
    row["order_relLinf"] = order_rel_linf[i]
    row["farField_order_relL2"] = order_ff_rel_l2[i]
    row["farField_order_relLinf"] = order_ff_rel_linf[i]

fieldnames = [
    "cellsPerWavelength",
    "h_over_lambda",
    "time",
    "relL2",
    "relLinf",
    "order_relL2",
    "order_relLinf",
    "absLinf",
    "farField_relL2",
    "farField_relLinf",
    "farField_order_relL2",
    "farField_order_relLinf",
    "farField_absLinf",
]
with metrics_path.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in metrics:
        writer.writerow(row)

def fmt_float(val, precision=3):
    try:
        v = float(val)
    except (TypeError, ValueError):
        return "--"
    if not np.isfinite(v):
        return "--"
    return f"{v:.{precision}e}"

def fmt_order(val):
    try:
        v = float(val)
    except (TypeError, ValueError):
        return "--"
    if not np.isfinite(v):
        return "--"
    return f"{v:.2f}"

tex_path = out_root / "nearField_convergence_table.tex"
with tex_path.open("w", encoding="utf-8") as f:
    f.write("\\begin{table}[htbp]\n")
    f.write("\\centering\n")
    f.write("\\begin{tabular}{ccccc}\n")
    f.write("\\hline\n")
    f.write("Cells per wavelength & $h/\\lambda$ & $E_2$ & $E_\\infty$ & Observed order \\\\\n")
    f.write("\\hline\n")
    for row in metrics:
        f.write(
            f"{row['cellsPerWavelength']} & "
            f"{float(row['h_over_lambda']):.4f} & "
            f"{fmt_float(row['relL2'])} & "
            f"{fmt_float(row['relLinf'])} & "
            f"{fmt_order(row['order_relL2'])} \\\\\n"
        )
    f.write("\\hline\n")
    f.write("\\end{tabular}\n")
    f.write("\\caption{Mesh-convergence study for the piston-radiation near-field pressure amplitude along the symmetry axis.}\n")
    f.write("\\label{tab:pistonConvergence}\n")
    f.write("\\end{table}\n")

ff_tex_path = out_root / "farField_convergence_table.tex"
with ff_tex_path.open("w", encoding="utf-8") as f:
    f.write("\\begin{table}[htbp]\n")
    f.write("\\centering\n")
    f.write("\\begin{tabular}{ccccc}\n")
    f.write("\\hline\n")
    f.write("Cells per wavelength & $h/\\lambda$ & $E_2^{ff}$ & $E_\\infty^{ff}$ & Observed order \\\\\n")
    f.write("\\hline\n")
    for row in metrics:
        f.write(
            f"{row['cellsPerWavelength']} & "
            f"{float(row['h_over_lambda']):.4f} & "
            f"{fmt_float(row['farField_relL2'])} & "
            f"{fmt_float(row['farField_relLinf'])} & "
            f"{fmt_order(row['farField_order_relL2'])} \\\\\n"
        )
    f.write("\\hline\n")
    f.write("\\end{tabular}\n")
    f.write("\\caption{Mesh-convergence study for the piston-radiation far-field directivity.}\n")
    f.write("\\label{tab:pistonFarFieldConvergence}\n")
    f.write("\\end{table}\n")

fig, ax = plt.subplots(figsize=(7.0, 4.6))
x = np.array([float(row["cellsPerWavelength"]) for row in metrics])
ax.loglog(x, [float(row["relL2"]) for row in metrics], "o-", label="near-field relL2")
ax.loglog(x, [float(row["relLinf"]) for row in metrics], "o--", label="near-field relLinf")
ax.loglog(x, [float(row["farField_relL2"]) for row in metrics], "s-", label="far-field relL2")
ax.loglog(x, [float(row["farField_relLinf"]) for row in metrics], "s--", label="far-field relLinf")
ax.set_xlabel("cells per wavelength")
ax.set_ylabel("relative error")
ax.grid(True, which="both", alpha=0.35)
ax.legend(loc="best")
ax.set_title("pistonRadiation mesh-convergence metrics")
fig.tight_layout()
fig.savefig(out_root / "metrics_meshConvergence.png", dpi=180)
PY

echo "Mesh convergence outputs written to: $outRoot"
echo "Summary: $outRoot/metrics_summary.csv"
