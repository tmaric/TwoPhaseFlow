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


def fit_inc_ref(x: np.ndarray, p: np.ndarray, k: float):
    M = np.column_stack((np.exp(-1j*k*x), np.exp(1j*k*x)))
    c, *_ = np.linalg.lstsq(M, p, rcond=None)
    return c[0], c[1]


def fit_forward(x: np.ndarray, p: np.ndarray, k: float, x0: float):
    v = np.exp(-1j*k*(x - x0))
    c = np.vdot(v, p) / np.vdot(v, v)
    return c


def fit_fwd_bwd(x: np.ndarray, p: np.ndarray, k: float, x0: float):
    M = np.column_stack((np.exp(-1j*k*(x - x0)), np.exp(1j*k*(x - x0))))
    c, *_ = np.linalg.lstsq(M, p, rcond=None)
    return c[0], c[1]


def analytical_pressure(
    x: np.ndarray,
    incident_amp: complex,
    x_int: float,
    k1: float,
    k2: float,
    r_th: float,
    t_th: float,
) -> np.ndarray:
    p = np.empty_like(x, dtype=complex)
    left = x <= x_int
    incident_at_interface = incident_amp*np.exp(-1j*k1*x_int)
    reflected_amp = r_th*incident_at_interface*np.exp(-1j*k1*x_int)
    transmitted_amp = t_th*incident_at_interface
    p[left] = incident_amp*(
        np.exp(-1j*k1*x[left])
    ) + reflected_amp*np.exp(1j*k1*x[left])
    p[~left] = transmitted_amp*np.exp(-1j*k2*(x[~left] - x_int))
    return p


def main():
    case_dir = pathlib.Path(__file__).resolve().parent
    p = parse_keyvals(case_dir / "caseParams.sh")

    f = p["DRIVE_F"]
    x_int = p["X_INTERFACE"]
    x_max = p["X_MAX"]
    pml_l = p["PML_L"]
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

    x, pc = probe_line(vtu, (0.005, y, z), (x_max - 0.005, y, z), 1200)

    left_mask = (x > 0.03) & (x < x_int - 0.015)
    right_mask = (x > x_int + 0.015) & (x < (x_max - pml_l - 0.015))
    if np.count_nonzero(left_mask) < 20 or np.count_nonzero(right_mask) < 20:
        raise RuntimeError("Not enough sampling points for fitting")

    A, B = fit_inc_ref(x[left_mask], pc[left_mask], k1)
    C, D = fit_fwd_bwd(x[right_mask], pc[right_mask], k2, x_int)

    incident_at_interface = A*np.exp(-1j*k1*x_int)
    reflected_at_interface = B*np.exp(1j*k1*x_int)
    R_num = reflected_at_interface/incident_at_interface
    T_num = C/incident_at_interface

    Z1 = rho1*c1
    Z2 = rho2*c2
    R_th = (Z2 - Z1)/(Z2 + Z1)
    T_th = (2*Z2)/(Z2 + Z1)
    pc_th = analytical_pressure(x, A, x_int, k1, k2, R_th, T_th)
    compare_mask = (x > 0.005) & (x < x_max - pml_l - 0.005)
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
        fobj.write(f"R_num = {R_num.real:.8e} + i{R_num.imag:.8e}\n")
        fobj.write(f"R_th  = {R_th:.8e}\n")
        fobj.write(f"|R_num|-|R_th| = {abs(abs(R_num)-abs(R_th)):.8e}\n")
        fobj.write(f"T_num = {T_num.real:.8e} + i{T_num.imag:.8e}\n")
        fobj.write(f"T_th  = {T_th:.8e}\n")
        fobj.write(f"|T_num|-|T_th| = {abs(abs(T_num)-abs(T_th)):.8e}\n")
        fobj.write(f"Right-going ratio D/C in medium2 = {abs(D/C):.8e}\n")
        fobj.write(f"Pre_relL2 = {pre_rel_l2:.8e}\n")
        fobj.write(f"Pim_relL2 = {pim_rel_l2:.8e}\n")
        fobj.write(f"|p|_relL2 = {pabs_rel_l2:.8e}\n")

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=120)
    labels = ["|R|", "|T|"]
    num = [abs(R_num), abs(T_num)]
    ana = [abs(R_th), abs(T_th)]
    ix = np.arange(len(labels))
    wbar = 0.35
    ax.bar(ix - 0.5*wbar, ana, wbar, label="Analytical")
    ax.bar(ix + 0.5*wbar, num, wbar, label="Simulation")
    ax.set_xticks(ix)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Coefficient magnitude")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "interface_RT_compare.png")

    fig_p, axes = plt.subplots(2, 1, figsize=(8, 6), dpi=140, sharex=True)
    axes[0].plot(x, pc_th.real, "k-", lw=1.8, label="Analytical Pre")
    axes[0].plot(x, pc.real, "r--", lw=1.3, label="Simulation Pre")
    axes[1].plot(x, pc_th.imag, "k-", lw=1.8, label="Analytical Pim")
    axes[1].plot(x, pc.imag, "b--", lw=1.3, label="Simulation Pim")
    for axi in axes:
        axi.axvline(x_int, color="0.35", lw=1.0, ls=":", label="interface")
        axi.axvline(x_max - pml_l, color="0.55", lw=1.0, ls="--", label="PML start")
        axi.set_ylabel("pressure [Pa]")
        axi.grid(True, alpha=0.3)
        handles, labels = axi.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        axi.legend(unique.values(), unique.keys(), loc="best")
    axes[1].set_xlabel("x [m]")
    fig_p.suptitle("Layered interface pressure field comparison")
    fig_p.tight_layout()
    fig_p.savefig(out_dir / "pressureField_Pre_Pim_compare.png")

    fig_abs, ax_abs = plt.subplots(figsize=(8, 4.5), dpi=140)
    ax_abs.plot(x, np.abs(pc_th), "k-", lw=1.8, label="Analytical |p|")
    ax_abs.plot(x, np.abs(pc), "g--", lw=1.3, label="Simulation |p|")
    ax_abs.axvline(x_int, color="0.35", lw=1.0, ls=":", label="interface")
    ax_abs.axvline(x_max - pml_l, color="0.55", lw=1.0, ls="--", label="PML start")
    ax_abs.set_xlabel("x [m]")
    ax_abs.set_ylabel("|p| [Pa]")
    ax_abs.grid(True, alpha=0.3)
    ax_abs.legend(loc="best")
    fig_abs.tight_layout()
    fig_abs.savefig(out_dir / "pressureField_abs_compare.png")

    print(f"Wrote: {out_dir/'metrics.txt'}")
    print(f"Wrote: {out_dir/'interface_RT_compare.png'}")
    print(f"Wrote: {pressure_csv}")
    print(f"Wrote: {out_dir/'pressureField_Pre_Pim_compare.png'}")
    print(f"Wrote: {out_dir/'pressureField_abs_compare.png'}")


if __name__ == "__main__":
    main()
