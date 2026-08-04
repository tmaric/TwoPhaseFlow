#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import re
import subprocess
import numpy as np

import matplotlib.pyplot as plt
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkFiltersCore import vtkProbeFilter
from vtkmodules.vtkFiltersCore import vtkCellDataToPointData
from vtkmodules.vtkFiltersSources import vtkLineSource
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader


def parse_keyvals(path: pathlib.Path) -> dict[str, float]:
    vals: dict[str, float] = {}
    pat = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([0-9eE+\-.]+)\s*$")
    for line in path.read_text(encoding="utf-8").splitlines():
        m = pat.match(line)
        if m:
            vals[m.group(1)] = float(m.group(2))
    return vals


def latest_time(case_dir: pathlib.Path) -> str:
    times = []
    for p in case_dir.iterdir():
        if p.is_dir():
            try:
                times.append((float(p.name), p.name))
            except ValueError:
                pass
    if not times:
        raise RuntimeError("No time directories found")
    return sorted(times)[-1][1]


def latest_vtu(case_dir: pathlib.Path, time_name: str) -> pathlib.Path:
    vtu = list((case_dir / "VTK").glob(f"*_{time_name}/internal.vtu"))
    if not vtu:
        raise RuntimeError("No internal.vtu produced")
    return vtu[0]


def probe_line(vtu: pathlib.Path, p1: tuple[float, float, float], p2: tuple[float, float, float], n: int):
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(vtu))
    reader.Update()

    c2p = vtkCellDataToPointData()
    c2p.SetInputData(reader.GetOutput())
    c2p.PassCellDataOn()
    c2p.Update()

    line = vtkLineSource()
    line.SetPoint1(*p1)
    line.SetPoint2(*p2)
    line.SetResolution(n - 1)
    line.Update()

    probe = vtkProbeFilter()
    probe.SetInputData(line.GetOutput())
    probe.SetSourceData(c2p.GetOutput())
    probe.Update()

    out = probe.GetOutput()
    pts = vtk_to_numpy(out.GetPoints().GetData())
    pre = vtk_to_numpy(out.GetPointData().GetArray("Pre"))
    pim = vtk_to_numpy(out.GetPointData().GetArray("Pim"))
    return pts[:, 0], pre + 1j*pim


def analytical_pressure(
    x: np.ndarray,
    rho1: float,
    c1: float,
    c2: float,
    u0: float,
    k1: float,
    k2: float,
    interface_x: float,
    pml_start_x: float,
    domain_end_x: float,
    sigma_max: float,
    pml_order: float,
    R: float,
) -> np.ndarray:
    """Closed-form layered-interface pressure including PML attenuation."""
    p = np.empty_like(x, dtype=complex)
    left = x < interface_x
    amplitude = rho1*c1*u0
    denominator = 1.0 - R*np.exp(2j*k1*interface_x)
    p[left] = amplitude*(
        np.exp(1j*k1*x[left])
        + R*np.exp(1j*k1*(2.0*interface_x - x[left]))
    )/denominator
    p[~left] = amplitude*(1.0 + R)*np.exp(
        1j*(k2*(x[~left] - interface_x) + k1*interface_x)
    )/denominator

    pml = x > pml_start_x
    pml_length = domain_end_x - pml_start_x
    normalized_distance = (x[pml] - pml_start_x)/pml_length
    attenuation = np.exp(
        sigma_max*pml_length/((pml_order + 1.0)*c2)
        * normalized_distance**(pml_order + 1.0)
    )
    p[pml] *= attenuation

    return p


def main():
    case_dir = pathlib.Path(__file__).resolve().parent
    p = parse_keyvals(case_dir / "studyProperties")

    f = p["DRIVE_F"]
    piston_u = p["PISTON_U"]
    x_int = p["X_INTERFACE"]
    x_max = p["X_MAX"]
    pml_l = p["PML_L"]
    pml_start = x_max - pml_l
    sigma_max = -p["SIGMA_MAX"]
    pml_order = p["PO"]
    y = 0.0
    z = 0.0

    rho1, c1 = p["RHOG"], p["CG"]
    rho2, c2 = p["RHOL"], p["CL"]
    w = 2*np.pi*f
    k1 = w/c1
    k2 = w/c2

    subprocess.run(["foamToVTK", "-case", str(case_dir), "-latestTime", "-noZero", "-fields", "(Pre Pim)"], check=True, stdout=subprocess.DEVNULL)

    t = latest_time(case_dir)
    vtu = latest_vtu(case_dir, t)

    x, pc = probe_line(vtu, (0.0, y, z), (x_max, y, z), 1200)

    Z1 = rho1*c1
    Z2 = rho2*c2
    R = (Z2 - Z1)/(Z2 + Z1)
    pc_th = analytical_pressure(
        x,
        rho1,
        c1,
        c2,
        piston_u,
        k1,
        k2,
        x_int,
        pml_start,
        x_max,
        sigma_max,
        pml_order,
        R,
    )
    # Use the complete gas--liquid--PML domain for the verification norm.
    compare_mask = np.isfinite(pc) & np.isfinite(pc_th)
    p_rel_l2 = np.linalg.norm((pc - pc_th)[compare_mask]) / max(
        np.linalg.norm(pc_th[compare_mask]), 1e-30
    )
    pre_rel_l2 = np.linalg.norm((pc.real - pc_th.real)[compare_mask]) / max(
        np.linalg.norm(pc_th.real[compare_mask]), 1e-30
    )
    pim_rel_l2 = np.linalg.norm((pc.imag - pc_th.imag)[compare_mask]) / max(
        np.linalg.norm(pc_th.imag[compare_mask]), 1e-30
    )
    pabs_rel_l2 = np.linalg.norm((np.abs(pc) - np.abs(pc_th))[compare_mask]) / max(
        np.linalg.norm(np.abs(pc_th)[compare_mask]), 1e-30
    )

    out_dir = case_dir / "postProcessing" / "interfaceValidation" / t
    out_dir.mkdir(parents=True, exist_ok=True)

    pressure_csv = out_dir / "pressureFieldComparison.csv"
    np.savetxt(
        pressure_csv,
        np.column_stack(
            [
                x,
                pc.real,
                pc.imag,
                np.abs(pc),
                pc_th.real,
                pc_th.imag,
                np.abs(pc_th),
            ]
        ),
        delimiter=",",
        header="x_m,Pre_sim,Pim_sim,p_abs_sim,Pre_analytic,Pim_analytic,p_abs_analytic",
        comments="",
    )

    with (out_dir / "metrics.txt").open("w", encoding="utf-8") as fobj:
        fobj.write("Layered interface validation\n")
        fobj.write("Analytical pressure field includes polynomial PML attenuation\n")
        fobj.write(f"rho1 = {rho1:.8e}\n")
        fobj.write(f"c1 = {c1:.8e}\n")
        fobj.write(f"rho2 = {rho2:.8e}\n")
        fobj.write(f"c2 = {c2:.8e}\n")
        fobj.write(f"u0 = {piston_u:.8e}\n")
        fobj.write(f"L1_interface = {x_int:.8e}\n")
        fobj.write(f"L2_PML_start = {pml_start:.8e}\n")
        fobj.write(f"L_domain_end = {x_max:.8e}\n")
        fobj.write(f"k1 = {k1:.8e}\n")
        fobj.write(f"k2 = {k2:.8e}\n")
        fobj.write(f"R = {R:.8e}\n")
        fobj.write(f"PML_length = {pml_l:.8e}\n")
        fobj.write(f"sigmaMax = {sigma_max:.8e}\n")
        fobj.write(f"PML_order = {pml_order:.8e}\n")
        fobj.write("errorDomain = whole domain including PML\n")
        fobj.write(f"P_relL2 = {p_rel_l2:.8e}\n")
        fobj.write(f"Pre_relL2 = {pre_rel_l2:.8e}\n")
        fobj.write(f"Pim_relL2 = {pim_rel_l2:.8e}\n")
        fobj.write(f"|p|_relL2 = {pabs_rel_l2:.8e}\n")

    fig_p, axes = plt.subplots(2, 1, figsize=(8, 6), dpi=180, sharex=True)
    components = [
        (pc_th.real, pc.real, r"$P_{\mathrm{re}}$", "#d55e00"),
        (pc_th.imag, pc.imag, r"$P_{\mathrm{im}}$", "#0072b2"),
    ]
    for row, (analytical, simulation, label, color) in enumerate(components):
        ax = axes[row]
        ax.plot(x, analytical, color="0.1", lw=2.0, label="Analytical")
        ax.plot(
            x,
            simulation,
            color=color,
            lw=1.2,
            ls="--",
            marker="o",
            markevery=45,
            ms=2.2,
            mfc="white",
            mec=color,
            mew=0.7,
            label="Numerical",
        )
        ax.axvspan(pml_start, x_max, color="0.92", zorder=0)
        ax.axvline(x_int, color="0.35", lw=0.9, ls=":", label="interface")
        ax.axvline(pml_start, color="0.55", lw=0.9, ls="--", label="PML start")
        ax.set_xlim(0.0, x_max)
        ax.set_ylabel(f"{label} [Pa]")
        ax.grid(True, alpha=0.25)
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        ax.legend(
            unique.values(),
            unique.keys(),
            loc="upper right",
            bbox_to_anchor=(0.985, 0.98),
            frameon=True,
            framealpha=0.9,
            facecolor="white",
            edgecolor="0.85",
            fancybox=False,
            ncol=1,
            fontsize="small",
            labelspacing=0.25,
            handlelength=1.8,
            borderaxespad=0.2,
        )

    axes[0].set_title("Layered interface pressure field comparison")
    axes[1].set_xlabel("x [m]")
    fig_p.tight_layout()
    fig_p.savefig(out_dir / "pressureField_Pre_Pim_compare.png")

    fig_abs, ax_abs = plt.subplots(figsize=(8, 4.5), dpi=180)
    ax_abs.plot(x, np.abs(pc_th), color="0.1", lw=2.0, label="Analytical |p|")
    ax_abs.plot(
        x,
        np.abs(pc),
        color="#009e73",
        lw=1.2,
        ls="--",
        marker="o",
        markevery=45,
        ms=2.2,
        mfc="white",
        mec="#009e73",
        mew=0.7,
        label="Numerical |p|",
    )
    ax_abs.axvspan(pml_start, x_max, color="0.92", zorder=0)
    ax_abs.axvline(x_int, color="0.35", lw=1.0, ls=":", label="interface")
    ax_abs.axvline(x_max - pml_l, color="0.55", lw=1.0, ls="--", label="PML start")
    ax_abs.set_xlim(0.0, x_max)
    ax_abs.set_title("Pressure amplitude comparison")
    ax_abs.set_xlabel("x [m]")
    ax_abs.set_ylabel("|p| [Pa]")
    ax_abs.grid(True, alpha=0.25)
    ax_abs.legend(
        loc="upper right",
        bbox_to_anchor=(0.985, 0.98),
        frameon=True,
        framealpha=0.9,
        facecolor="white",
        edgecolor="0.85",
        fancybox=False,
        ncol=1,
        fontsize="small",
        labelspacing=0.25,
        handlelength=1.8,
        borderaxespad=0.2,
    )
    fig_abs.tight_layout()
    fig_abs.savefig(out_dir / "pressureField_abs_compare.png")

    print(f"Wrote: {out_dir/'metrics.txt'}")
    print(f"Wrote: {pressure_csv}")
    print(f"Wrote: {out_dir/'pressureField_Pre_Pim_compare.png'}")
    print(f"Wrote: {out_dir/'pressureField_abs_compare.png'}")


if __name__ == "__main__":
    main()
