#!/usr/bin/env python3
"""
Compare pistonRadiation simulation with analytical baffled-piston on-axis pressure.

Analytical reference (Rayleigh integral on axis):
    p(z) = rho*c*u0 * (exp(-i*k*z) - exp(-i*k*sqrt(z^2 + a^2)))
"""

from __future__ import annotations

import argparse
import glob
import math
import pathlib
import re
import subprocess
import sys

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkFiltersCore import vtkProbeFilter
from vtkmodules.vtkFiltersSources import vtkLineSource
from vtkmodules.vtkFiltersCore import vtkCellDataToPointData
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader

try:
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover
    raise SystemExit(f"matplotlib is required: {exc}")


def parse_case_params(path: pathlib.Path) -> dict[str, float]:
    wanted = {"DRIVE_F", "PISTON_U", "CG", "RHOG", "PML_RMIN"}
    out: dict[str, float] = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            m = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([0-9eE+.\-]+)\s*$", line)
            if not m:
                continue
            key, val = m.group(1), m.group(2)
            if key in wanted:
                out[key] = float(val)
    missing = wanted - set(out)
    if missing:
        raise RuntimeError(f"Missing in {path}: {sorted(missing)}")
    return out


def parse_piston_radius(path: pathlib.Path) -> float:
    txt = path.read_text(encoding="utf-8")
    m = re.search(r"radius\s+([0-9eE+.\-]+)\s*;\s*\}\s*\}\s*\)", txt)
    if m:
        return float(m.group(1))
    m = re.search(r"\bradius\s+([0-9eE+.\-]+)\s*;", txt)
    if not m:
        raise RuntimeError(f"Could not parse piston radius from {path}")
    return float(m.group(1))


def latest_time(case_dir: pathlib.Path) -> str:
    times = []
    for p in case_dir.iterdir():
        if not p.is_dir():
            continue
        try:
            val = float(p.name)
        except ValueError:
            continue
        times.append((val, p.name))
    if not times:
        raise RuntimeError("No time directories found in case")
    times.sort()
    return times[-1][1]


def ensure_case_ready(case_dir: pathlib.Path) -> None:
    points = case_dir / "constant" / "polyMesh" / "points"
    if not points.exists():
        raise RuntimeError(
            f"Mesh not found at {points}. Run ./Allrun first."
        )
    _ = latest_time(case_dir)


def run_foam_to_vtk(case_dir: pathlib.Path) -> None:
    cmd = [
        "foamToVTK",
        "-case",
        str(case_dir),
        "-latestTime",
        "-noZero",
        "-fields",
        "(Pre Pim)",
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def latest_internal_vtu(case_dir: pathlib.Path, time_name: str) -> pathlib.Path:
    vtk_dir = case_dir / "VTK"
    candidates = sorted(glob.glob(str(vtk_dir / f"*_{time_name}" / "internal.vtu")))
    if not candidates:
        raise RuntimeError(f"No VTK internal file found for time {time_name} under {vtk_dir}")
    return pathlib.Path(candidates[0])


def extract_axis_profile_from_vtu(
    vtu_path: pathlib.Path, z_max: float, n_points: int, a: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(vtu_path))
    reader.Update()
    ug = reader.GetOutput()

    c2p = vtkCellDataToPointData()
    c2p.SetInputData(ug)
    c2p.PassCellDataOn()
    c2p.Update()
    src = c2p.GetOutput()

    line = vtkLineSource()
    line.SetPoint1(0.0, 1e-6, 0.0)
    line.SetPoint2(0.0, z_max, 0.0)
    line.SetResolution(max(1, n_points - 1))
    line.Update()

    probe = vtkProbeFilter()
    probe.SetInputData(line.GetOutput())
    probe.SetSourceData(src)
    probe.Update()
    out = probe.GetOutput()

    pts = vtk_to_numpy(out.GetPoints().GetData())
    y = pts[:, 1]

    pre_arr = out.GetPointData().GetArray("Pre")
    pim_arr = out.GetPointData().GetArray("Pim")
    if pre_arr is None or pim_arr is None:
        raise RuntimeError(f"Probed Pre/Pim not found in {vtu_path}")
    pre = vtk_to_numpy(pre_arr)
    pim = vtk_to_numpy(pim_arr)

    valid = np.isfinite(pre) & np.isfinite(pim) & (y > 0)
    if np.count_nonzero(valid) < 4:
        raise RuntimeError("Insufficient valid probed axis points.")

    return y[valid], pre[valid], pim[valid]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", default=".", help="OpenFOAM case directory")
    parser.add_argument("--n-points", type=int, default=400, help="Sampling points on axis")
    parser.add_argument("--z-max", type=float, default=None, help="Axis max distance [m]")
    parser.add_argument("--a", type=float, default=None, help="Piston radius [m] override")
    args = parser.parse_args()

    case_dir = pathlib.Path(args.case).resolve()
    ensure_case_ready(case_dir)
    params = parse_case_params(case_dir / "caseParams.sh")

    f = params["DRIVE_F"]
    u0 = params["PISTON_U"]
    c = params["CG"]
    rho = params["RHOG"]
    a = args.a if args.a is not None else parse_piston_radius(case_dir / "system" / "topoSetDict")
    z_max = args.z_max if args.z_max is not None else 0.9 * params["PML_RMIN"]

    t_name = latest_time(case_dir)
    run_foam_to_vtk(case_dir)
    vtu = latest_internal_vtu(case_dir, t_name)
    z, pre, pim = extract_axis_profile_from_vtu(vtu, z_max, args.n_points, a)

    k = 2.0 * math.pi * f / c
    p0 = rho * c * u0
    r = np.sqrt(z * z + a * a)
    p_analytic = p0 * np.abs(np.exp(-1j * k * z) - np.exp(-1j * k * r))
    p_sim = np.sqrt(pre * pre + pim * pim)

    rayleigh = 0.5 * k * a * a
    x = z / rayleigh
    y_sim = p_sim / p0
    y_ana = p_analytic / p0

    err_l2 = float(np.linalg.norm(y_sim - y_ana) / max(np.linalg.norm(y_ana), 1e-30))
    err_linf = float(np.max(np.abs(y_sim - y_ana)))

    out_dir = case_dir / "postProcessing" / "analyticalCompare" / t_name
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "onAxisComparison.csv"
    hdr = "z_m,z_over_rayleigh,p_sim_pa,p_analytic_pa,p_sim_over_p0,p_analytic_over_p0"
    np.savetxt(
        csv_path,
        np.column_stack([z, x, p_sim, p_analytic, y_sim, y_ana]),
        delimiter=",",
        header=hdr,
        comments="",
    )

    fig = plt.figure(figsize=(7.2, 4.6))
    ax = fig.add_subplot(1, 1, 1)
    # Match COMSOL-style axis display: log scale z/R0 from 1e-3 to 5e0.
    x_plot = np.clip(x, 1e-12, None)
    ax.plot(x_plot, y_ana, "-", lw=2.0, label="Analytical (Rayleigh on-axis)")
    ax.plot(x_plot, y_sim, "--", lw=1.8, label="Simulation (|Pre+iPim|)")
    ax.set_xlabel("z / R0")
    ax.set_ylabel("|p| / (rho*c*u0)")
    ax.set_xscale("log")
    ax.set_xlim(5e-4, 1)
    ax.grid(True, alpha=0.35)
    ax.legend(loc="best")
    ax.set_title(f"pistonRadiation on-axis comparison, f={f:.0f} Hz")
    fig.tight_layout()
    png_path = out_dir / "onAxisComparison.png"
    fig.savefig(png_path, dpi=160)

    report = out_dir / "metrics.txt"
    report.write_text(
        "\n".join(
            [
                f"time = {t_name}",
                f"f_Hz = {f}",
                f"a_m = {a}",
                f"rayleigh_m = {rayleigh}",
                f"relL2 = {err_l2:.6e}",
                f"absLinf = {err_linf:.6e}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Wrote: {csv_path}")
    print(f"Wrote: {png_path}")
    print(f"Wrote: {report}")
    print(f"Metrics: relL2={err_l2:.6e}, absLinf={err_linf:.6e}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as err:
        print(f"Error: {err}", file=sys.stderr)
        raise SystemExit(1)
