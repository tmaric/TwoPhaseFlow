#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import pathlib
import re
import subprocess

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkFiltersCore import vtkCellCenters
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader


def properties(path: pathlib.Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def latest_time(case_dir: pathlib.Path) -> str:
    times = []
    for item in case_dir.iterdir():
        if item.is_dir():
            try:
                times.append((float(item.name), item.name))
            except ValueError:
                pass
    return max(times)[1]


def mesh_quality(log_path: pathlib.Path) -> tuple[float, float]:
    text = log_path.read_text(encoding="utf-8", errors="replace")
    match = re.search(
        r"Mesh non-orthogonality Max:\s*([-+0-9.eE]+)\s+average:\s*([-+0-9.eE]+)",
        text,
    )
    if not match:
        return float("nan"), float("nan")
    return float(match.group(1)), float(match.group(2))


def rel_l2(error: np.ndarray, reference: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sqrt(np.sum(weights*error**2) / np.sum(weights*reference**2)))


def main() -> None:
    case_dir = pathlib.Path(__file__).resolve().parent
    p = properties(case_dir / "studyProperties")
    subprocess.run(
        ["foamToVTK", "-case", str(case_dir), "-latestTime", "-noZero", "-fields", "(Pre Pim Ure Uim V)"],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    time_name = latest_time(case_dir)
    candidates = list((case_dir / "VTK").glob(f"*_{time_name}/internal.vtu"))
    if not candidates:
        raise RuntimeError("foamToVTK did not create internal.vtu")

    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(candidates[0]))
    reader.Update()
    grid = reader.GetOutput()
    centres_filter = vtkCellCenters()
    centres_filter.SetInputData(grid)
    centres_filter.VertexCellsOn()
    centres_filter.Update()
    centres = vtk_to_numpy(centres_filter.GetOutput().GetPoints().GetData())

    data = grid.GetCellData()
    pre = vtk_to_numpy(data.GetArray("Pre"))
    pim = vtk_to_numpy(data.GetArray("Pim"))
    ure = vtk_to_numpy(data.GetArray("Ure"))
    uim = vtk_to_numpy(data.GetArray("Uim"))
    volumes_array = data.GetArray("V")
    volumes = vtk_to_numpy(volumes_array) if volumes_array else np.ones_like(pre)

    f = float(p["DRIVE_F"])
    rho = float(p["RHO0"])
    c = float(p["C0"])
    amp = float(p["P0"])
    phase = float(p["PHASE"])
    k = 2*math.pi*f/c
    if p["BOUNDARY_MODE"] == "dirichlet":
        kvec = np.array([k*math.cos(math.pi/6), k*math.sin(math.pi/6), 0.0])
    else:
        kvec = np.array([k, 0.0, 0.0])
    theta = phase - centres @ kvec
    pre_exact = amp*np.cos(theta)
    pim_exact = amp*np.sin(theta)
    ure_exact = -(amp/(2*math.pi*f*rho))*np.cos(theta)[:, None]*kvec
    uim_exact = -(amp/(2*math.pi*f*rho))*np.sin(theta)[:, None]*kvec

    complex_error = np.sqrt((pre-pre_exact)**2 + (pim-pim_exact)**2)
    complex_reference = np.sqrt(pre_exact**2 + pim_exact**2)
    velocity_error = np.sqrt(np.sum((ure-ure_exact)**2 + (uim-uim_exact)**2, axis=1))
    velocity_reference = np.sqrt(np.sum(ure_exact**2 + uim_exact**2, axis=1))

    lx = float(p["LX"])
    ly = float(p["LY"])
    h = 1.0/float(p["CELLS_PER_WAVELENGTH"])
    distance = np.minimum.reduce(
        [centres[:, 0], lx-centres[:, 0], centres[:, 1], ly-centres[:, 1]]
    )
    interior = distance > 2.5*h
    boundary = ~interior
    max_nonorth, avg_nonorth = mesh_quality(case_dir / "log.checkMesh")

    metrics = {
        "cellsPerWavelength": float(p["CELLS_PER_WAVELENGTH"]),
        "cells": float(grid.GetNumberOfCells()),
        "hOverLambda": h,
        "maxNonOrthogonality": max_nonorth,
        "avgNonOrthogonality": avg_nonorth,
        "pressureRelL2": rel_l2(complex_error, complex_reference, volumes),
        "pressureRelLinf": float(np.max(complex_error)/np.max(complex_reference)),
        "velocityRelL2": rel_l2(velocity_error, velocity_reference, volumes),
        "velocityInteriorRelL2": rel_l2(
            velocity_error[interior], velocity_reference[interior], volumes[interior]
        ) if np.any(interior) else float("nan"),
        "velocityBoundaryRelL2": rel_l2(
            velocity_error[boundary], velocity_reference[boundary], volumes[boundary]
        ),
    }

    out = case_dir / "postProcessing" / "homogeneousBaseline" / time_name
    out.mkdir(parents=True, exist_ok=True)
    with (out / "metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics))
        writer.writeheader()
        writer.writerow(metrics)
    with (out / "metrics.txt").open("w", encoding="utf-8") as stream:
        for key, value in metrics.items():
            stream.write(f"{key} = {value:.12e}\n")
    print(out / "metrics.txt")


if __name__ == "__main__":
    main()
