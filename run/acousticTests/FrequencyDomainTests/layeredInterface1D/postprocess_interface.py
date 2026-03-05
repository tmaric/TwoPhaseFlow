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

    R_num = B/A
    T_num = C/A

    Z1 = rho1*c1
    Z2 = rho2*c2
    R_th = (Z2 - Z1)/(Z2 + Z1)
    T_th = (2*Z2)/(Z2 + Z1)

    out_dir = case_dir / "postProcessing" / "interfaceValidation" / t
    out_dir.mkdir(parents=True, exist_ok=True)

    with (out_dir / "metrics.txt").open("w", encoding="utf-8") as fobj:
        fobj.write("Layered interface validation\n")
        fobj.write(f"R_num = {R_num.real:.8e} + i{R_num.imag:.8e}\n")
        fobj.write(f"R_th  = {R_th:.8e}\n")
        fobj.write(f"|R_num|-|R_th| = {abs(abs(R_num)-abs(R_th)):.8e}\n")
        fobj.write(f"T_num = {T_num.real:.8e} + i{T_num.imag:.8e}\n")
        fobj.write(f"T_th  = {T_th:.8e}\n")
        fobj.write(f"|T_num|-|T_th| = {abs(abs(T_num)-abs(T_th)):.8e}\n")
        fobj.write(f"Right-going ratio D/C in medium2 = {abs(D/C):.8e}\n")

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

    print(f"Wrote: {out_dir/'metrics.txt'}")
    print(f"Wrote: {out_dir/'interface_RT_compare.png'}")


if __name__ == "__main__":
    main()
