#!/usr/bin/env python3
"""
Compare pistonRadiation with the exact baffled-piston pressure field.

The reference is the Rayleigh integral of a uniform piston in a rigid baffle,
evaluated as plane wave minus edge wave (piston_reference.py, validated by
validate_reference.py).  Cell-centre values are compared with the reference at
the same cell centres, without interpolation:

- field: complex pressure in the physical region (outside the PML),
  volume-weighted relative L2 error;
- axis: complex pressure in the cell row next to the symmetry axis.

The pressure level on an arc of radius r_arc is written for plotting; these
values are interpolated by VTK probing and do not enter the error metrics.
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
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData, vtkVertex
from vtkmodules.vtkFiltersCore import vtkCellDataToPointData, vtkProbeFilter
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader

from piston_reference import piston_pressure, piston_pressure_on_axis

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover
    raise SystemExit(f"matplotlib is required: {exc}")


def parse_case_params(path: pathlib.Path) -> dict[str, float]:
    wanted = {
        "DRIVE_F",
        "PISTON_U",
        "PISTON_RADIUS",
        "CG",
        "RHOG",
        "PML_RMIN",
    }
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


def mesh_cell_count(case_dir: pathlib.Path) -> int:
    with (case_dir / "constant" / "polyMesh" / "owner").open("rb") as f:
        header = f.read(4096).decode("latin-1")
    m = re.search(r"nCells:\s*(\d+)", header)
    if not m:
        raise RuntimeError("nCells not found in the polyMesh/owner header")
    return int(m.group(1))


def write_vtk_with_cell_geometry(case_dir: pathlib.Path) -> None:
    """Write cell centres and volumes, then convert the latest time to VTK."""
    for func in ("writeCellCentres", "writeCellVolumes"):
        subprocess.run(
            ["postProcess", "-case", str(case_dir), "-func", func, "-latestTime"],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
    subprocess.run(
        [
            "foamToVTK", "-case", str(case_dir), "-latestTime", "-noZero",
            "-fields", "(Pre Pim C V)",
        ],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )


def latest_internal_vtu(case_dir: pathlib.Path, time_name: str) -> pathlib.Path:
    vtk_dir = case_dir / "VTK"
    candidates = sorted(glob.glob(str(vtk_dir / f"*_{time_name}" / "internal.vtu")))
    if not candidates:
        raise RuntimeError(f"No VTK internal file found for time {time_name} under {vtk_dir}")
    return pathlib.Path(candidates[0])


def read_unstructured_grid(vtu_path: pathlib.Path):
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(vtu_path))
    reader.Update()
    return reader.GetOutput()


def cell_arrays(ug, names: tuple[str, ...]) -> list[np.ndarray]:
    out = []
    for name in names:
        arr = ug.GetCellData().GetArray(name)
        if arr is None:
            raise RuntimeError(f"Cell array {name} not found in the VTK output")
        out.append(vtk_to_numpy(arr))
    return out


def probe_points_pre_pim(ug, pts_xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    c2p = vtkCellDataToPointData()
    c2p.SetInputData(ug)
    c2p.PassCellDataOn()
    c2p.Update()

    vtk_pts = vtkPoints()
    verts = vtkCellArray()
    for p in pts_xyz:
        pid = vtk_pts.InsertNextPoint(float(p[0]), float(p[1]), float(p[2]))
        v = vtkVertex()
        v.GetPointIds().SetId(0, pid)
        verts.InsertNextCell(v)

    pdata = vtkPolyData()
    pdata.SetPoints(vtk_pts)
    pdata.SetVerts(verts)

    probe = vtkProbeFilter()
    probe.SetInputData(pdata)
    probe.SetSourceData(c2p.GetOutput())
    probe.Update()
    out = probe.GetOutput()

    pre_arr = out.GetPointData().GetArray("Pre")
    pim_arr = out.GetPointData().GetArray("Pim")
    if pre_arr is None or pim_arr is None:
        raise RuntimeError("Probed Pre/Pim not found")
    return vtk_to_numpy(pre_arr), vtk_to_numpy(pim_arr)


def spl_db(p_abs: np.ndarray, p_ref: float = 20e-6) -> np.ndarray:
    # Sound pressure level of a harmonic amplitude |P|: 20 log10(|P|/(sqrt(2) p_ref))
    return 20.0 * np.log10(np.maximum(p_abs, 1e-30) / (math.sqrt(2.0) * p_ref))


def rel_l2(num: np.ndarray, ref: np.ndarray, weight: np.ndarray | None = None) -> float:
    w = np.ones(num.shape) if weight is None else weight
    return float(np.sqrt(np.sum(w * np.abs(num - ref) ** 2) / np.sum(w * np.abs(ref) ** 2)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", default=".", help="OpenFOAM case directory")
    parser.add_argument("--n-theta", type=int, default=361, help="Samples on the plotting arc")
    parser.add_argument(
        "--r-arc",
        type=float,
        default=None,
        help="Plotting-arc radius [m] (default: 0.9*PML_RMIN)",
    )
    parser.add_argument(
        "--n-rim",
        type=int,
        default=1024,
        help="Gauss-Legendre nodes of the rim integral in the reference",
    )
    parser.add_argument(
        "--output-name",
        default="analyticalCompare",
        help="Postprocessing output directory name",
    )
    args = parser.parse_args()

    case_dir = pathlib.Path(args.case).resolve()
    ensure_case_ready(case_dir)
    params = parse_case_params(case_dir / "caseParams.sh")

    f = params["DRIVE_F"]
    c = params["CG"]
    a = params["PISTON_RADIUS"]
    extent = params["PML_RMIN"]
    r_arc = args.r_arc if args.r_arc is not None else 0.9 * extent
    k = 2.0 * math.pi * f / c
    p0 = params["RHOG"] * c * params["PISTON_U"]
    rayleigh = 0.5 * k * a * a

    t_name = latest_time(case_dir)
    write_vtk_with_cell_geometry(case_dir)
    ug = read_unstructured_grid(latest_internal_vtu(case_dir, t_name))
    pre, pim, C, V = cell_arrays(ug, ("Pre", "Pim", "C", "V"))
    n_cells = mesh_cell_count(case_dir)
    if len(V) != n_cells:
        raise RuntimeError(f"VTK has {len(V)} cells, the mesh {n_cells} (decomposed polyhedra?)")

    # Axis y, radial x, wedge z; the PML starts at PML_RMIN in r and y.
    r = np.hypot(C[:, 0], C[:, 2])
    z = C[:, 1]
    phys = (r <= extent) & (z <= extent)
    p_h = pre[phys] + 1j * pim[phys]
    p_ex = piston_pressure(r[phys], z[phys], k, a, p0, n_phi=args.n_rim)
    r_p, z_p, V_p = r[phys], z[phys], V[phys]

    err_field = rel_l2(p_h, p_ex, V_p)
    err_field_mag = rel_l2(np.abs(p_h), np.abs(p_ex), V_p)

    axis = r_p <= r_p.min() * (1.0 + 1e-9)
    order = np.argsort(z_p[axis])
    z_ax, r_ax = z_p[axis][order], r_p[axis][order]
    ph_ax, pe_ax = p_h[axis][order], p_ex[axis][order]
    err_axis = rel_l2(ph_ax, pe_ax)
    err_axis_inf = float(np.max(np.abs(ph_ax - pe_ax)) / np.max(np.abs(pe_ax)))
    err_axis_mag = rel_l2(np.abs(ph_ax), np.abs(pe_ax))

    out_dir = case_dir / "postProcessing" / args.output_name / t_name
    out_dir.mkdir(parents=True, exist_ok=True)

    axis_csv = out_dir / "onAxisComparison.csv"
    np.savetxt(
        axis_csv,
        np.column_stack([
            z_ax, r_ax, z_ax / rayleigh,
            ph_ax.real, ph_ax.imag, pe_ax.real, pe_ax.imag,
            np.abs(ph_ax) / p0, np.abs(pe_ax) / p0,
        ]),
        delimiter=",",
        header="z_m,r_m,z_over_rayleigh,Pre_sim_pa,Pim_sim_pa,Pre_exact_pa,Pim_exact_pa,"
               "p_sim_over_p0,p_exact_over_p0",
        comments="",
    )

    theta = np.linspace(0.0, 0.5 * math.pi, max(3, args.n_theta))
    x_arc, z_arc = r_arc * np.sin(theta), r_arc * np.cos(theta)
    pre_arc, pim_arc = probe_points_pre_pim(
        ug, np.column_stack([x_arc, z_arc, np.zeros_like(theta)])
    )
    p_arc = pre_arc + 1j * pim_arc
    p_arc_ex = piston_pressure(x_arc, z_arc, k, a, p0, n_phi=args.n_rim)
    arc_csv = out_dir / "arcComparison.csv"
    np.savetxt(
        arc_csv,
        np.column_stack([
            np.degrees(theta), pre_arc, pim_arc, p_arc_ex.real, p_arc_ex.imag,
            spl_db(np.abs(p_arc)), spl_db(np.abs(p_arc_ex)), np.full_like(theta, r_arc),
        ]),
        delimiter=",",
        header="theta_deg,Pre_sim_pa,Pim_sim_pa,Pre_exact_pa,Pim_exact_pa,"
               "SPL_sim_dB,SPL_exact_dB,r_arc_m",
        comments="",
    )

    field_csv = out_dir / "fieldComparison.csv.gz"
    np.savetxt(
        field_csv,
        np.column_stack([r_p, z_p, V_p, p_h.real, p_h.imag, p_ex.real, p_ex.imag]),
        delimiter=",",
        header="r_m,z_m,V_m3,Pre_sim_pa,Pim_sim_pa,Pre_exact_pa,Pim_exact_pa",
        comments="",
    )

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    z_line = np.geomspace(1e-4, extent, 2000)
    ax.plot(z_line / rayleigh, np.abs(piston_pressure_on_axis(z_line, k, a, p0)) / p0,
            "-", lw=2.0, label="Exact (axis)")
    ax.plot(z_ax / rayleigh, np.abs(ph_ax) / p0, "--", lw=1.6, label="Numerical (cell row at the axis)")
    ax.set_xscale("log")
    ax.set_xlabel("z / R0")
    ax.set_ylabel("|p| / (rho*c*u0)")
    ax.grid(True, alpha=0.35)
    ax.legend(loc="best")
    ax.set_title(f"pistonRadiation on-axis comparison, f={f:.0f} Hz")
    fig.tight_layout()
    axis_png = out_dir / "onAxisComparison.png"
    fig.savefig(axis_png, dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(np.degrees(theta), spl_db(np.abs(p_arc_ex)), "-", lw=2.0, label="Exact")
    ax.plot(np.degrees(theta), spl_db(np.abs(p_arc)), "--", lw=1.8, label="Numerical")
    ax.set_xlabel(r"$\theta$ [$^\circ$] from the symmetry axis")
    ax.set_ylabel(r"$L_p$ [dB re $20\,\mu$Pa]")
    ax.set_xlim(0, 90)
    ax.grid(True, alpha=0.35)
    ax.legend(loc="best")
    ax.set_title(f"Pressure level on the arc r = {r_arc:g} m")
    fig.tight_layout()
    arc_png = out_dir / "arcComparison.png"
    fig.savefig(arc_png, dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.0, 5.2))
    err = np.log10(np.maximum(np.abs(p_h - p_ex) / p0, 1e-12))
    tpc = ax.tripcolor(r_p, z_p, err, shading="flat", cmap="viridis")
    fig.colorbar(tpc, ax=ax, label=r"$\log_{10}(|p_h-p|/(\rho c u_0))$")
    ax.plot([0, a], [0, 0], color="red", lw=3)
    ax.set_xlabel("r [m]")
    ax.set_ylabel("z [m]")
    ax.set_aspect("equal")
    ax.set_title("Pointwise error in the physical region")
    fig.tight_layout()
    err_png = out_dir / "fieldError.png"
    fig.savefig(err_png, dpi=160)
    plt.close(fig)

    report = out_dir / "metrics.txt"
    report.write_text(
        "\n".join(
            [
                f"time = {t_name}",
                f"f_Hz = {f}",
                f"a_m = {a}",
                f"rho_c_u0_pa = {p0}",
                f"rayleigh_m = {rayleigh}",
                f"physical_extent_m = {extent}",
                f"r_arc_m = {r_arc}",
                f"rim_nodes = {args.n_rim}",
                f"field_cells = {int(phys.sum())}",
                f"axis_cells = {int(axis.sum())}",
                f"field_relL2 = {err_field:.6e}",
                f"field_magnitude_relL2 = {err_field_mag:.6e}",
                f"axis_relL2 = {err_axis:.6e}",
                f"axis_relLinf = {err_axis_inf:.6e}",
                f"axis_magnitude_relL2 = {err_axis_mag:.6e}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    for path in (axis_csv, arc_csv, field_csv, axis_png, arc_png, err_png, report):
        print(f"Wrote: {path}")
    print(
        "Metrics: "
        f"field(relL2={err_field:.6e}) "
        f"axis(relL2={err_axis:.6e}, relLinf={err_axis_inf:.6e})"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as err:
        print(f"Error: {err}", file=sys.stderr)
        raise SystemExit(1)
