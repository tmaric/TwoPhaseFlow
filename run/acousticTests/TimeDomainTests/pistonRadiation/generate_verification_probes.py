#!/usr/bin/env python3
import json
import math
import re
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
SYSTEM_DIR = CASE_DIR / "system"
TRANSPORT = CASE_DIR / "constant" / "transportProperties"
P_FIELD = CASE_DIR / "0" / "p"


def read_scalar(path: Path, key: str) -> float:
    text = path.read_text(encoding="utf-8", errors="ignore")
    # match: key   value; or key value;
    m = re.search(rf"^\s*{re.escape(key)}\s+([-+0-9.eE]+)\s*;", text, re.MULTILINE)
    if not m:
        raise ValueError(f"Missing {key} in {path}")
    return float(m.group(1))


def main() -> int:
    a = read_scalar(TRANSPORT, "pistonRadius")
    R0 = read_scalar(TRANSPORT, "domainRadius")

    # Sampling plan (small set)
    n_axis = 7
    n_theta = 9
    y_start = 0.5 * a
    y_end = 4.0 * a
    r_ring = 0.5 * R0

    # Axis points along +y (x=z=0)
    axis_pts = []
    for i in range(n_axis):
        y = y_start + (y_end - y_start) * i / (n_axis - 1)
        axis_pts.append((0.0, y, 0.0))

    # Ring points at fixed r from origin, varying polar angle from axis
    # theta in [0, pi/2] for the wedge (x>=0, y>=0, z=0)
    ring_pts = []
    for i in range(n_theta):
        theta = (math.pi / 2) * i / (n_theta - 1)
        x = r_ring * math.sin(theta)
        y = r_ring * math.cos(theta)
        ring_pts.append((x, y, 0.0))

    all_pts = axis_pts + ring_pts

    # Write include file for controlDict (list items only, no surrounding parens)
    loc_path = SYSTEM_DIR / "verificationProbesLocations"
    with loc_path.open("w", encoding="utf-8") as f:
        for x, y, z in all_pts:
            f.write(f"    ({x:.8g} {y:.8g} {z:.8g})\n")

    meta = {
        "axis_indices": list(range(0, len(axis_pts))),
        "ring_indices": list(range(len(axis_pts), len(all_pts))),
        "r_ring": r_ring,
        "n_axis": n_axis,
        "n_theta": n_theta,
        "y_start": y_start,
        "y_end": y_end,
    }
    meta_path = SYSTEM_DIR / "verificationProbesMeta.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"Wrote {loc_path}")
    print(f"Wrote {meta_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
