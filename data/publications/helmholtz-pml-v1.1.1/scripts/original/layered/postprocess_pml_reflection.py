#!/usr/bin/env python3
from __future__ import annotations

import csv
import pathlib
import subprocess

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkFiltersCore import vtkCellCenters
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader


def properties(path: pathlib.Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        try:
            values[key] = float(value)
        except ValueError:
            pass
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
    p = properties(case_dir / "studyProperties")
    subprocess.run(
        ["foamToVTK", "-case", str(case_dir), "-latestTime", "-noZero", "-fields", "(Pre Pim)"],
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

    omega = 2*np.pi*p["DRIVE_F"]
    k = omega/p["CG"]
    wavelength = p["CG"]/p["DRIVE_F"]
    pml_start = p["X_MAX"]-p["PML_L"]
    fit = (x > p["X_MIN"]+wavelength) & (x < pml_start-0.5*wavelength)
    basis = np.column_stack([np.exp(1j*k*x[fit]), np.exp(-1j*k*x[fit])])
    incident, reflected = np.linalg.lstsq(basis, pressure[fit], rcond=None)[0]
    fitted = basis @ np.array([incident, reflected])
    fit_residual = np.linalg.norm(pressure[fit]-fitted)/np.linalg.norm(pressure[fit])

    pml = x >= pml_start
    start_index = np.argmin(np.abs(x-pml_start))
    end_index = np.argmax(x)
    attenuation = abs(pressure[end_index])/max(abs(pressure[start_index]), 1e-30)
    metrics = {
        "NX": p["NX"],
        "cellsPerWavelength": p["NX"]*wavelength/(p["X_MAX"]-p["X_MIN"]),
        "sigmaOverOmega": p["SIGMA_MAX"]/omega,
        "pmlThicknessOverLambda": p["PML_L"]/wavelength,
        "pmlOrder": p["PO"],
        "reflectionMagnitude": abs(reflected/incident),
        "fitResidual": fit_residual,
        "outerToInnerAmplitude": attenuation,
    }
    out = case_dir / "postProcessing" / "pmlReflection" / time_name
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

