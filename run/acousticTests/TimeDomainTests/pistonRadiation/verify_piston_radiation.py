#!/usr/bin/env python3
import argparse
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
from mpmath import besselj


def j1(x: float) -> float:
    return float(besselj(1, x))

CASE_DIR = Path(__file__).resolve().parent
SYSTEM_DIR = CASE_DIR / "system"
TRANSPORT = CASE_DIR / "constant" / "transportProperties"
P_FIELD = CASE_DIR / "0" / "p"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Verify piston radiation vs analytical.")
    p.add_argument("--plots", action="store_true", help="Generate PNG plots via gnuplot")
    return p.parse_args()


def read_scalar(path: Path, key: str) -> float:
    text = path.read_text(encoding="utf-8", errors="ignore")
    m = re.search(
        rf"^\s*{re.escape(key)}\s+(?:\[[^\]]+\]\s+)?([-+0-9.eE]+)\s*;",
        text,
        re.MULTILINE,
    )
    if not m:
        raise ValueError(f"Missing {key} in {path}")
    return float(m.group(1))


def latest_probe_file() -> Path:
    probes_dir = CASE_DIR / "postProcessing" / "verificationProbes"
    if not probes_dir.exists():
        raise FileNotFoundError(f"Missing {probes_dir}")
    subdirs = [d for d in probes_dir.iterdir() if d.is_dir()]
    if not subdirs:
        raise FileNotFoundError(f"No probe output in {probes_dir}")
    subdirs.sort(key=lambda p: float(p.name))
    for d in reversed(subdirs):
        fpath = d / "p"
        if fpath.exists():
            return fpath
    raise FileNotFoundError("No probe file 'p' found")


def read_probe_timeseries(path: Path):
    times = []
    values = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                t = float(parts[0])
                row = [float(x) for x in parts[1:]]
            except ValueError:
                continue
            times.append(t)
            values.append(row)
    return np.array(times), np.array(values)


def amplitude_over_window(times, values, t_end, window):
    t_start = t_end - window
    mask = (times >= t_start) & (times <= t_end)
    if not np.any(mask):
        raise ValueError("No samples in window")
    v = values[mask, :]
    vmax = np.max(v, axis=0)
    vmin = np.min(v, axis=0)
    return 0.5 * (vmax - vmin)


def near_field_on_axis_amp(r, a, rho0, c0, U0, k):
    # From piston_radiation.pdf: p(r,0) = 2*rho0*c0*U0 * sin( (k r /2) * (sqrt(1+(a/r)^2) - 1) )
    term = (k * r / 2.0) * (math.sqrt(1.0 + (a / r) ** 2) - 1.0)
    return 2.0 * rho0 * c0 * U0 * abs(math.sin(term))


def far_field_amp(r, theta, a, rho0, U0, omega, k):
    x = k * a * math.sin(theta)
    if abs(x) < 1e-12:
        dir_factor = 1.0
    else:
        dir_factor = 2.0 * j1(x) / x
    return (omega * rho0 * a * a * U0) / (2.0 * r) * abs(dir_factor)


def main() -> int:
    args = parse_args()
    # Read constants
    a = read_scalar(TRANSPORT, "pistonRadius")
    R0 = read_scalar(TRANSPORT, "domainRadius")
    c0 = read_scalar(TRANSPORT, "soundSpeed")
    rho0 = read_scalar(TRANSPORT, "rho0")
    f = read_scalar(TRANSPORT, "f")
    U0 = read_scalar(P_FIELD, "pistonVelocity")

    omega = 2.0 * math.pi * f
    k = omega / c0
    period = 1.0 / f
    window = 5.0 * period

    meta_path = SYSTEM_DIR / "verificationProbesMeta.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))

    probe_path = latest_probe_file()
    times, values = read_probe_timeseries(probe_path)
    if times.size == 0:
        raise ValueError("Empty probe file")

    amp = amplitude_over_window(times, values, times[-1], window)

    # Load probe locations
    loc_path = SYSTEM_DIR / "verificationProbesLocations"
    locs = []
    for line in loc_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line == "(" or line == ")":
            continue
        line = line.strip("()")
        parts = line.split()
        if len(parts) != 3:
            continue
        locs.append(tuple(float(x) for x in parts))

    # Axis comparison
    axis_rows = []
    for idx in meta["axis_indices"]:
        x, y, z = locs[idx]
        r = math.sqrt(x * x + y * y + z * z)
        p_sim = amp[idx]
        p_ana = near_field_on_axis_amp(r, a, rho0, c0, U0, k)
        rel = abs(p_sim - p_ana) / p_ana if p_ana != 0 else 0.0
        axis_rows.append((y, p_sim, p_ana, rel))

    # Ring comparison (far-field)
    ring_rows = []
    for idx in meta["ring_indices"]:
        x, y, z = locs[idx]
        r = math.sqrt(x * x + y * y + z * z)
        theta = math.atan2(math.sqrt(x * x + z * z), y)
        p_sim = amp[idx]
        p_ana = far_field_amp(r, theta, a, rho0, U0, omega, k)
        rel = abs(p_sim - p_ana) / p_ana if p_ana != 0 else 0.0
        ring_rows.append((theta, p_sim, p_ana, rel))

    # Write CSVs
    out_axis = CASE_DIR / "verification_axis.csv"
    out_ring = CASE_DIR / "verification_ring.csv"
    with out_axis.open("w", encoding="utf-8") as f:
        f.write("y,p_sim,p_ana,rel_err\n")
        for y, p_sim, p_ana, rel in axis_rows:
            f.write(f"{y},{p_sim},{p_ana},{rel}\n")
    with out_ring.open("w", encoding="utf-8") as f:
        f.write("theta_rad,theta_deg,p_sim,p_ana,rel_err\n")
        for theta, p_sim, p_ana, rel in ring_rows:
            f.write(f"{theta},{math.degrees(theta)},{p_sim},{p_ana},{rel}\n")

    # Optional plots (off by default)
    if args.plots:
        gnuplot = shutil.which("gnuplot")
        if gnuplot:
            axis_png = CASE_DIR / "verification_axis.png"
            ring_png = CASE_DIR / "verification_ring.png"
            subprocess.run(
                [
                    gnuplot,
                    "-e",
                    (
                        "set datafile separator ','; "
                        "set terminal png size 900,600; "
                        f"set output '{axis_png}'; "
                        "set xlabel 'y [m]'; "
                        "set ylabel 'pressure amplitude [Pa]'; "
                        "set key left top; "
                        f"plot '{out_axis}' using 1:2 with linespoints title 'sim', "
                        f"'' using 1:3 with linespoints title 'analytical';"
                    ),
                ],
                check=False,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            subprocess.run(
                [
                    gnuplot,
                    "-e",
                    (
                        "set datafile separator ','; "
                        "set terminal png size 900,600; "
                        f"set output '{ring_png}'; "
                        "set xlabel 'polar angle [deg]'; "
                        "set ylabel 'pressure amplitude [Pa]'; "
                        "set key left top; "
                        f"plot '{out_ring}' using 2:3 with linespoints title 'sim', "
                        f"'' using 2:4 with linespoints title 'analytical';"
                    ),
                ],
                check=False,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

    print(f"Wrote {out_axis}")
    print(f"Wrote {out_ring}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
