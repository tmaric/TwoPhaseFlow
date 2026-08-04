#!/usr/bin/env python3
from __future__ import annotations

import csv
import pathlib

import matplotlib.pyplot as plt


def scientific_tex(value: float) -> str:
    mantissa, exponent = f"{value:.3e}".split("e")
    return "$" + mantissa + r"\times 10^{" + str(int(exponent)) + "}$"


def main() -> None:
    case_dir = pathlib.Path(__file__).resolve().parent
    out_dir = case_dir / "postProcessing" / "pmlMeshSensitivity"
    csv_path = out_dir / "metrics_summary.csv"

    rows = []
    with csv_path.open(newline="", encoding="utf-8") as fobj:
        reader = csv.DictReader(fobj)
        for row in reader:
            rows.append(
                {
                    "sigmaMax": float(row["sigmaMax"]),
                    "NX": int(float(row["NX"])),
                    "dx_m": float(row["dx_m"]),
                    "P_relL2": float(row["P_relL2"]),
                    "Pre_relL2": float(row["Pre_relL2"]),
                    "Pim_relL2": float(row["Pim_relL2"]),
                }
            )
    rows.sort(key=lambda item: item["NX"])

    with csv_path.open("w", newline="", encoding="utf-8") as fobj:
        writer = csv.DictWriter(
            fobj,
            fieldnames=[
                "sigmaMax",
                "NX",
                "dx_m",
                "P_relL2",
                "Pre_relL2",
                "Pim_relL2",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "sigmaMax": f"{row['sigmaMax']:.0f}",
                    "NX": f"{row['NX']}",
                    "dx_m": f"{row['dx_m']:.8e}",
                    "P_relL2": f"{row['P_relL2']:.8e}",
                    "Pre_relL2": f"{row['Pre_relL2']:.8e}",
                    "Pim_relL2": f"{row['Pim_relL2']:.8e}",
                }
            )

    fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=180)
    ax.loglog(
        [row["NX"] for row in rows],
        [row["P_relL2"] for row in rows],
        "o-",
        lw=2.0,
        label=r"$E_P$",
    )
    ax.loglog(
        [row["NX"] for row in rows],
        [row["Pre_relL2"] for row in rows],
        "s--",
        lw=1.3,
        label=r"$E_{P_{\mathrm{re}}}$",
    )
    ax.loglog(
        [row["NX"] for row in rows],
        [row["Pim_relL2"] for row in rows],
        "^:",
        lw=1.5,
        label=r"$E_{P_{\mathrm{im}}}$",
    )
    ax.set_xlabel(r"number of cells $N_x$")
    ax.set_ylabel(r"relative $L_2$ error")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_dir / "pml_mesh_sensitivity.png")

    sigma = rows[0]["sigmaMax"] if rows else 0.0
    sigma_text = f"{abs(sigma)/1.0e5:.1f}\\times10^5"

    with (out_dir / "pml_mesh_sensitivity_table.tex").open("w", encoding="utf-8") as fobj:
        fobj.write("\\begin{table}[H]\n")
        fobj.write("    \\centering\n")
        fobj.write(
            "    \\caption{Directional tensor PML mesh sensitivity for the "
            "one-dimensional layered-interface case with fixed damping "
            f"$\\sigma_{{\\max}}={sigma_text}\\,\\mathrm{{s^{{-1}}}}$.}}\n"
        )
        fobj.write("    \\label{tab:pmlMeshSensitivity}\n")
        fobj.write("    \\begin{tabular}{rrrrr}\n")
        fobj.write("        \\toprule\n")
        fobj.write(
            "        $N_x$ & $\\Delta x$ [$\\mathrm{m}$] & "
            "$E_P$ & $E_{P_{\\mathrm{re}}}$ & $E_{P_{\\mathrm{im}}}$ \\\\\n"
        )
        fobj.write("        \\midrule\n")
        for row in rows:
            fobj.write(
                "        "
                f"{row['NX']} & "
                f"{scientific_tex(row['dx_m'])} & "
                f"{scientific_tex(row['P_relL2'])} & "
                f"{scientific_tex(row['Pre_relL2'])} & "
                f"{scientific_tex(row['Pim_relL2'])} \\\\\n"
            )
        fobj.write("        \\bottomrule\n")
        fobj.write("    \\end{tabular}\n")
        fobj.write("\\end{table}\n")


if __name__ == "__main__":
    main()
