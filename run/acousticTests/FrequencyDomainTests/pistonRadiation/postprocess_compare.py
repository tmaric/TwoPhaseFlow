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
from scipy.special import j1
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkFiltersCore import vtkProbeFilter
from vtkmodules.vtkFiltersSources import vtkLineSource
from vtkmodules.vtkFiltersCore import vtkCellDataToPointData
from vtkmodules.vtkFiltersGeneral import vtkGradientFilter
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData, vtkVertex
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader

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


def probe_points_pre_pim(vtu_path: pathlib.Path, pts_xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(vtu_path))
    reader.Update()
    ug = reader.GetOutput()

    c2p = vtkCellDataToPointData()
    c2p.SetInputData(ug)
    c2p.PassCellDataOn()
    c2p.Update()
    src = c2p.GetOutput()

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
    probe.SetSourceData(src)
    probe.Update()
    out = probe.GetOutput()

    pre_arr = out.GetPointData().GetArray("Pre")
    pim_arr = out.GetPointData().GetArray("Pim")
    if pre_arr is None or pim_arr is None:
        raise RuntimeError(f"Probed Pre/Pim not found in {vtu_path}")
    pre = vtk_to_numpy(pre_arr)
    pim = vtk_to_numpy(pim_arr)
    return pre, pim


def probe_points_pre_pim_and_grad(
    vtu_path: pathlib.Path, pts_xyz: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(vtu_path))
    reader.Update()
    ug = reader.GetOutput()

    c2p = vtkCellDataToPointData()
    c2p.SetInputData(ug)
    c2p.PassCellDataOn()
    c2p.Update()
    src = c2p.GetOutput()

    g_pre = vtkGradientFilter()
    g_pre.SetInputData(src)
    g_pre.SetInputScalars(0, "Pre")
    g_pre.SetResultArrayName("gradPre")
    g_pre.Update()

    g_pim = vtkGradientFilter()
    g_pim.SetInputData(g_pre.GetOutput())
    g_pim.SetInputScalars(0, "Pim")
    g_pim.SetResultArrayName("gradPim")
    g_pim.Update()

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
    probe.SetSourceData(g_pim.GetOutput())
    probe.Update()
    out = probe.GetOutput()

    pre_arr = out.GetPointData().GetArray("Pre")
    pim_arr = out.GetPointData().GetArray("Pim")
    grad_pre_arr = out.GetPointData().GetArray("gradPre")
    grad_pim_arr = out.GetPointData().GetArray("gradPim")
    if pre_arr is None or pim_arr is None or grad_pre_arr is None or grad_pim_arr is None:
        raise RuntimeError(f"Probed Pre/Pim or gradients not found in {vtu_path}")

    pre = vtk_to_numpy(pre_arr)
    pim = vtk_to_numpy(pim_arr)
    grad_pre = vtk_to_numpy(grad_pre_arr)
    grad_pim = vtk_to_numpy(grad_pim_arr)
    return pre, pim, grad_pre, grad_pim


def far_field_directivity_analytic(theta: np.ndarray, k: float, a: float) -> np.ndarray:
    x = k * a * np.sin(theta)
    out = np.ones_like(x)
    small = np.abs(x) < 1e-10
    out[~small] = 2.0 * j1(x[~small]) / x[~small]
    out[small] = 1.0
    return np.abs(out)


def safe_spl_db(p_abs: np.ndarray, p_ref: float = 20e-6) -> np.ndarray:
    # Match COMSOL expression: 10*log10(0.5*|p|^2/p_ref^2)
    ratio = 0.5 * (np.maximum(p_abs, 1e-30) ** 2) / (p_ref**2)
    return 10.0 * np.log10(np.maximum(ratio, 1e-300))


def latest_openair_vtp(case_dir: pathlib.Path, time_name: str) -> pathlib.Path:
    vtk_dir = case_dir / "VTK"
    candidates = sorted(glob.glob(str(vtk_dir / f"*_{time_name}" / "boundary" / "openAir.vtp")))
    if not candidates:
        raise RuntimeError(f"No openAir.vtp found for time {time_name} under {vtk_dir}")
    return pathlib.Path(candidates[0])


def load_openair_boundary(vtp_path: pathlib.Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r = vtkXMLPolyDataReader()
    r.SetFileName(str(vtp_path))
    r.Update()
    d = r.GetOutput()
    pts = vtk_to_numpy(d.GetPoints().GetData())

    pre_arr = d.GetCellData().GetArray("Pre")
    pim_arr = d.GetCellData().GetArray("Pim")
    if pre_arr is None or pim_arr is None:
        raise RuntimeError(f"Cell Pre/Pim not found in {vtp_path}")
    pre = vtk_to_numpy(pre_arr)
    pim = vtk_to_numpy(pim_arr)
    p = pre + 1j * pim

    n_cells = d.GetNumberOfCells()
    centers = np.zeros((n_cells, 3))
    normals = np.zeros((n_cells, 3))
    areas = np.zeros(n_cells)

    for ci in range(n_cells):
        cell = d.GetCell(ci)
        nids = cell.GetNumberOfPoints()
        if nids < 3:
            continue
        pids = [cell.GetPointId(j) for j in range(nids)]
        cpts = pts[pids, :]
        centers[ci, :] = np.mean(cpts, axis=0)

        # Polygon area-vector by fan triangulation
        p0 = cpts[0]
        area_vec = np.zeros(3)
        for j in range(1, nids - 1):
            v1 = cpts[j] - p0
            v2 = cpts[j + 1] - p0
            area_vec += 0.5 * np.cross(v1, v2)

        area = float(np.linalg.norm(area_vec))
        if area > 0.0:
            n = area_vec / area
            # Ensure outward orientation (away from origin for this case)
            if np.dot(n, centers[ci, :]) < 0.0:
                n = -n
            normals[ci, :] = n
        areas[ci] = area

    valid = areas > 0.0
    return centers[valid], normals[valid], areas[valid], p[valid]


def build_inner_pml_source_from_vtu(
    vtu_path: pathlib.Path, r_src: float, n_theta: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # Build a virtual spherical source contour in axisymmetric plane:
    # x = r_src*sin(theta), y = r_src*cos(theta), theta in [0, pi/2]
    theta = np.linspace(0.0, 0.5 * math.pi, max(8, n_theta))
    x = r_src * np.sin(theta)
    y = r_src * np.cos(theta)
    pts = np.column_stack([x, y, np.zeros_like(theta)])

    pre, pim, grad_pre, grad_pim = probe_points_pre_pim_and_grad(vtu_path, pts)
    p = pre + 1j * pim

    # Outward normal for spherical contour
    n = pts / max(r_src, 1e-30)
    grad_p = grad_pre + 1j * grad_pim
    dpdn = np.einsum("ij,ij->i", grad_p, n)

    # Axisymmetric surface area represented by each theta sample
    dtheta = np.gradient(theta)
    dS = 2.0 * math.pi * (r_src**2) * np.sin(theta) * dtheta

    valid = (
        np.isfinite(p)
        & np.isfinite(dpdn)
        & np.isfinite(dS)
        & (dS > 0.0)
    )
    return pts[valid], n[valid], dS[valid], p[valid], dpdn[valid]


def revolve_axisymmetric_source(
    src_xyz_meridian: np.ndarray,
    src_area_meridian: np.ndarray,
    src_p_meridian: np.ndarray,
    src_dpdn_meridian: np.ndarray,
    n_phi: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # Build full 3D source quadrature by revolving meridian samples around y-axis.
    # Meridian coordinates are (x=r*sin(theta), y=r*cos(theta), z=0).
    if n_phi < 4:
        raise RuntimeError("n_phi must be >= 4 for full-phi integration")

    x_m = src_xyz_meridian[:, 0]
    y_m = src_xyz_meridian[:, 1]
    # Recover sphere radius from geometry.
    r_src = float(np.median(np.sqrt(np.sum(src_xyz_meridian**2, axis=1))))
    if r_src <= 0.0:
        raise RuntimeError("Invalid source radius while revolving axisymmetric contour")

    dphi = 2.0 * math.pi / float(n_phi)
    phi = (np.arange(n_phi, dtype=float) + 0.5) * dphi
    cphi = np.cos(phi)
    sphi = np.sin(phi)

    # Expand meridian data to (n_theta, n_phi) then flatten.
    x = np.outer(x_m, cphi)
    z = np.outer(x_m, sphi)
    y = np.outer(y_m, np.ones_like(phi))

    # Outward normal on sphere is radial unit vector.
    nx = x / r_src
    ny = y / r_src
    nz = z / r_src

    area = np.outer(src_area_meridian, np.full(n_phi, dphi / (2.0 * math.pi)))
    p = np.outer(src_p_meridian, np.ones(n_phi, dtype=complex))
    dpdn = np.outer(src_dpdn_meridian, np.ones(n_phi, dtype=complex))

    src_xyz = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
    src_n = np.column_stack([nx.ravel(), ny.ravel(), nz.ravel()])
    src_area = area.ravel()
    src_p = p.ravel()
    src_dpdn = dpdn.ravel()

    valid = (
        np.isfinite(src_area)
        & (src_area > 0.0)
        & np.isfinite(src_p)
        & np.isfinite(src_dpdn)
    )
    return src_xyz[valid], src_n[valid], src_area[valid], src_p[valid], src_dpdn[valid]


def kirchhoff_exterior_pressure(
    obs_xyz: np.ndarray,
    src_xyz: np.ndarray,
    src_n: np.ndarray,
    src_area: np.ndarray,
    src_p: np.ndarray,
    src_dpdn: np.ndarray,
    k: float,
) -> np.ndarray:
    # Full Kirchhoff integral:
    # p(P) = \int_S [ p*dG/dn - G*dp/dn ] dS
    out = np.zeros(obs_xyz.shape[0], dtype=np.complex128)
    for i in range(obs_xyz.shape[0]):
        rv = obs_xyz[i, :] - src_xyz
        R = np.linalg.norm(rv, axis=1)
        R = np.maximum(R, 1e-12)
        rhat = rv / R[:, None]
        G = np.exp(1j * k * R) / (4.0 * math.pi * R)
        dGdn = (1.0 / R - 1j * k) * G * np.einsum("ij,ij->i", rhat, src_n)
        out[i] = np.sum((src_p * dGdn - G * src_dpdn) * src_area)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", default=".", help="OpenFOAM case directory")
    parser.add_argument("--n-points", type=int, default=400, help="Sampling points on axis")
    parser.add_argument("--n-theta", type=int, default=361, help="Angular samples for far-field")
    parser.add_argument(
        "--n-phi",
        type=int,
        default=180,
        help="Azimuthal samples for full-phi source integration",
    )
    parser.add_argument("--z-max", type=float, default=None, help="Axis max distance [m]")
    parser.add_argument(
        "--r-ff1",
        type=float,
        default=None,
        help="(unused, kept for compatibility)",
    )
    parser.add_argument(
        "--r-ff2",
        type=float,
        default=None,
        help="(unused, kept for compatibility)",
    )
    parser.add_argument(
        "--r-far",
        type=float,
        default=None,
        help="Evaluation radius [m] for exterior far-field pattern (default: 20*PML_RMIN)",
    )
    parser.add_argument(
        "--r-src",
        type=float,
        default=None,
        help="Source contour radius [m] for boundary integral (default: 0.65*PML_RMIN)",
    )
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
    r_far = args.r_far if args.r_far is not None else 20.0 * params["PML_RMIN"]
    # Keep the reconstruction contour away from the PML coefficient transition,
    # where cell-to-point interpolation and numerical gradients are less robust.
    r_src = args.r_src if args.r_src is not None else 0.65*params["PML_RMIN"]

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

    abs_err = np.abs(y_sim - y_ana)
    err_l2 = float(np.linalg.norm(abs_err) / max(np.linalg.norm(y_ana), 1e-30))
    err_linf = float(np.max(abs_err) / max(float(np.max(np.abs(y_ana))), 1e-30))
    err_linf_abs_norm = float(np.max(abs_err))

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

    # Far-field/directivity comparison from exterior-field boundary integral
    # (Kirchhoff + Sommerfeld closure) evaluated on a large-radius arc.
    theta = np.linspace(0.0, 0.5 * math.pi, max(3, args.n_theta))
    x_far = r_far * np.sin(theta)
    y_far = r_far * np.cos(theta)
    pts_far = np.column_stack([x_far, y_far, np.zeros_like(theta)])

    src_xyz_m, _src_n_m, src_area_m, src_p_m, src_dpdn_m = build_inner_pml_source_from_vtu(
        vtu, r_src, max(32, args.n_theta)
    )
    src_xyz, src_n, src_area, src_p, src_dpdn = revolve_axisymmetric_source(
        src_xyz_m, src_area_m, src_p_m, src_dpdn_m, args.n_phi
    )
    p_far = kirchhoff_exterior_pressure(
        pts_far, src_xyz, src_n, src_area, src_p, src_dpdn, k
    )

    p_far_abs = np.abs(p_far)
    d_sim = p_far_abs.copy()
    d_sim /= max(float(d_sim[0]), 1e-30)  # normalize on-axis
    d_ana = far_field_directivity_analytic(theta, k, a)

    # Far-field analytical pressure magnitude from baffled piston theory.
    p_on_axis_ana = rho * c * u0 * k * a * a / (2.0 * r_far)
    p_ana_far_abs = p_on_axis_ana * d_ana
    spl_sim = safe_spl_db(p_far_abs)
    spl_ana = safe_spl_db(p_ana_far_abs)

    ff_abs_err = np.abs(d_sim - d_ana)
    err_ff_l2 = float(np.linalg.norm(ff_abs_err) / max(np.linalg.norm(d_ana), 1e-30))
    err_ff_linf = float(np.max(ff_abs_err) / max(float(np.max(np.abs(d_ana))), 1e-30))
    err_ff_linf_abs_norm = float(np.max(ff_abs_err))

    ff_csv = out_dir / "farFieldPatternComparison.csv"
    ff_hdr = (
        "theta_deg,SPL_sim_dB,SPL_analytic_dB,directivity_sim,directivity_analytic,"
        "p_far_abs_sim_pa,p_far_abs_analytic_pa,p_far_real,p_far_imag,r_far_m"
    )
    np.savetxt(
        ff_csv,
        np.column_stack(
            [
                np.degrees(theta),
                spl_sim,
                spl_ana,
                d_sim,
                d_ana,
                p_far_abs,
                p_ana_far_abs,
                np.real(p_far),
                np.imag(p_far),
                np.full_like(theta, r_far),
            ]
        ),
        delimiter=",",
        header=ff_hdr,
        comments="",
    )

    fig_ff = plt.figure(figsize=(7.2, 4.6))
    ax_ff = fig_ff.add_subplot(1, 1, 1)
    ax_ff.plot(np.degrees(theta), spl_ana, "-", lw=2.0, label="Analytical far-field SPL")
    ax_ff.plot(np.degrees(theta), spl_sim, "--", lw=1.8, label="Simulation SPL")
    ax_ff.set_xlabel("theta [deg] (from axis)")
    ax_ff.set_ylabel("SPL [dB re 20uPa]")
    ax_ff.set_xlim(0, 90)
    ax_ff.grid(True, alpha=0.35)
    ax_ff.legend(loc="best")
    ax_ff.set_title(f"pistonRadiation far-field SPL, f={f:.0f} Hz")
    fig_ff.tight_layout()
    ff_png = out_dir / "farFieldPatternComparison.png"
    fig_ff.savefig(ff_png, dpi=160)

    # COMSOL-style radiation pattern view (polar directivity).
    # Mirror the 0..90 deg axisymmetric result to show +/-90 deg lobe.
    th_m = np.concatenate((-theta[:0:-1], theta))
    spl_sim_m = np.concatenate((spl_sim[:0:-1], spl_sim))
    spl_ana_m = np.concatenate((spl_ana[:0:-1], spl_ana))

    fig_pol = plt.figure(figsize=(6.0, 6.0))
    ax_pol = fig_pol.add_subplot(1, 1, 1, projection="polar")
    ax_pol.set_theta_zero_location("N")
    ax_pol.set_theta_direction(-1)
    # Plot relative SPL in polar view (offset by each curve max for readability).
    rel_spl_ana = spl_ana_m - np.max(spl_ana_m)
    rel_spl_sim = spl_sim_m - np.max(spl_sim_m)
    ax_pol.plot(th_m, rel_spl_ana, "-", lw=2.0, label="Analytical")
    ax_pol.plot(th_m, rel_spl_sim, "--", lw=1.8, label="Simulation")
    ax_pol.set_thetamin(-90)
    ax_pol.set_thetamax(90)
    ax_pol.set_rlim(-40.0, 0.0)
    ax_pol.set_rlabel_position(105)
    ax_pol.grid(True, alpha=0.35)
    ax_pol.set_title("Far-field radiation pattern (relative SPL, dB)", pad=4)
    ax_pol.legend(loc="upper right", bbox_to_anchor=(1.03, 0.92))
    pol_png = out_dir / "farFieldPatternPolar.png"
    fig_pol.tight_layout()
    fig_pol.savefig(pol_png, dpi=160)

    report = out_dir / "metrics.txt"
    report.write_text(
        "\n".join(
            [
                f"time = {t_name}",
                f"f_Hz = {f}",
                f"a_m = {a}",
                f"rayleigh_m = {rayleigh}",
                f"r_src_m = {r_src}",
                f"r_far_m = {r_far}",
                f"relL2 = {err_l2:.6e}",
                f"relLinf = {err_linf:.6e}",
                f"absLinf = {err_linf_abs_norm:.6e}",
                f"farField_relL2 = {err_ff_l2:.6e}",
                f"farField_relLinf = {err_ff_linf:.6e}",
                f"farField_absLinf = {err_ff_linf_abs_norm:.6e}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Wrote: {csv_path}")
    print(f"Wrote: {png_path}")
    print(f"Wrote: {ff_csv}")
    print(f"Wrote: {ff_png}")
    print(f"Wrote: {pol_png}")
    print(f"Wrote: {report}")
    print(
        "Metrics: "
        f"onAxis(relL2={err_l2:.6e}, relLinf={err_linf:.6e}) "
        f"farFieldDir(relL2={err_ff_l2:.6e}, relLinf={err_ff_linf:.6e})"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as err:
        print(f"Error: {err}", file=sys.stderr)
        raise SystemExit(1)
