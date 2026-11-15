#!/usr/bin/env python3
from __future__ import annotations

import csv
import pathlib
import re
import subprocess

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkFiltersCore import vtkCellCenters
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader


def parameters(path: pathlib.Path) -> dict[str, float]:
    values: dict[str, float] = {}
    pattern = re.compile(r"^([A-Z][A-Z0-9_]*)=\$?\{?[^:}]*:-?([-+0-9.eE]+)\}?$")
    simple = re.compile(r"^([A-Z][A-Z0-9_]*)=([-+0-9.eE]+)$")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = simple.match(line.strip()) or pattern.match(line.strip())
        if match:
            values[match.group(1)] = float(match.group(2))
    return values


def latest_time(case_dir: pathlib.Path) -> str:
    times = []
    for path in case_dir.iterdir():
        if path.is_dir():
            try:
                times.append((float(path.name), path.name))
            except ValueError:
                pass
    return max(times)[1]


def main() -> None:
    case_dir = pathlib.Path(__file__).resolve().parent
    p = parameters(case_dir / "studyProperties")
    subprocess.run(
        ["foamToVTK", "-case", str(case_dir), "-latestTime", "-noZero", "-fields", "(Pre Pim V)"],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    time_name = latest_time(case_dir)
    vtu = list((case_dir / "VTK").glob(f"*_{time_name}/internal.vtu"))[0]
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(vtu))
    reader.Update()
    grid = reader.GetOutput()
    centres_filter = vtkCellCenters()
    centres_filter.SetInputData(grid)
    centres_filter.Update()
    x = vtk_to_numpy(centres_filter.GetOutput().GetPoints().GetData())[:, 0]
    data = grid.GetCellData()
    pressure = vtk_to_numpy(data.GetArray("Pre")) + 1j*vtk_to_numpy(data.GetArray("Pim"))
    volume_array = data.GetArray("V")
    volume = vtk_to_numpy(volume_array) if volume_array else np.ones_like(x)

    omega = 2*np.pi*p["DRIVE_F"]
    k1, k2 = omega/p["CG"], omega/p["CL"]
    xi = p["X_INTERFACE"]
    dx = (p["X_MAX"]-p["X_MIN"])/p["NX"]
    z1, z2 = p["RHOG"]*p["CG"], p["RHOL"]*p["CL"]
    r_exact = (z2-z1)/(z2+z1)
    t_exact = 1+r_exact

    left = (x > p["X_MIN"]+4*dx) & (x < xi-4*dx)
    right = (x > xi+4*dx) & (x < p["X_MAX"]-4*dx)
    left_basis = np.column_stack(
        [np.exp(1j*k1*(x[left]-xi)), np.exp(-1j*k1*(x[left]-xi))]
    )
    incident, reflected = np.linalg.lstsq(left_basis, pressure[left], rcond=None)[0]
    transmitted = np.linalg.lstsq(
        np.exp(1j*k2*(x[right]-xi))[:, None], pressure[right], rcond=None
    )[0][0]
    r_fit = reflected/incident
    t_fit = transmitted/incident

    exact = np.empty_like(pressure)
    exact[x < xi] = (
        np.exp(1j*k1*(x[x < xi]-xi))
        + r_exact*np.exp(-1j*k1*(x[x < xi]-xi))
    )
    exact[x >= xi] = t_exact*np.exp(1j*k2*(x[x >= xi]-xi))
    error = np.abs(pressure-exact)
    rel_l2 = np.sqrt(np.sum(volume*error**2)/np.sum(volume*np.abs(exact)**2))

    metrics = {
        "NX": p["NX"],
        "dx": dx,
        "interfaceX": xi,
        "pressureRelL2": rel_l2,
        "RReal": r_fit.real,
        "RImag": r_fit.imag,
        "RAbsError": abs(r_fit-r_exact),
        "TReal": t_fit.real,
        "TImag": t_fit.imag,
        "TAbsError": abs(t_fit-t_exact),
    }
    out = case_dir / "postProcessing" / "interfaceConvergence" / time_name
    out.mkdir(parents=True, exist_ok=True)
    with (out / "metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics))
        writer.writeheader()
        writer.writerow(metrics)
    with (out / "metrics.txt").open("w", encoding="utf-8") as stream:
        for key, value in metrics.items():
            stream.write(f"{key} = {value:.12e}\n")


if __name__ == "__main__":
    main()
