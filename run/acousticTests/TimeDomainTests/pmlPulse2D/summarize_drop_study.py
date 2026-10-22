#!/usr/bin/env python3
"""Tabulate the water-drop study E written by runStudyDrop.

Accuracy: E2 of the PML run against the enlarged-domain run on the same
cells over the physical domain, time-integrated up to T_obs, with energy
weights V beta for pressure and V rho for velocity; for all cells and for
the gas, liquid and mixed cells separately. Long time: physical energy and
whole-domain maxima as in study B.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

T_OBS = 4e-4


def probes(case, field):
    rows = [line.replace("(", " ").replace(")", " ").split() for line in
            (case / "postProcessing/cellProbes/0" / field).read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    data = np.array(rows, dtype=float)
    return data[:, 0], data[:, 1:]


def media(case):
    s = json.loads((case / "caseSetup.json").read_text())
    n, h, half, a = s["n"], s["h"], s["domainHalfWidth"], s["halfWidth"]
    text = (case / "0" / "alpha.water").read_text()
    body = text.split("internalField", 1)[1].split(";", 1)[0]
    values = body.split("(", 1)[1].rsplit(")", 1)[0].split()
    alpha = np.array(values, dtype=float).reshape(n, n)
    inside = np.abs(-half + (np.arange(n) + 0.5) * h) < a
    alpha = alpha[np.ix_(inside, inside)].ravel()
    rho = alpha * s["rhoL"] + (1 - alpha) * s["rhoG"]
    beta = alpha / (s["rhoL"] * s["cL"]**2) + (1 - alpha) / (s["rhoG"] * s["cG"]**2)
    return s, alpha, rho, beta


def e2(field, ref, weight, mask):
    w = weight[None, mask]
    if field.ndim == 3:
        w = w[..., None]
        diff, base = field[:, mask] - ref[:, mask], ref[:, mask]
    else:
        diff, base = field[:, mask] - ref[:, mask], ref[:, mask]
    return math.sqrt(np.sum(w * diff**2) / np.sum(w * base**2))


def accuracy(root):
    rows, prev = [], None
    for case in sorted((root / "accuracy").glob("pml_h*"),
                       key=lambda c: -float(c.name[5:])):
        ref = case.parent / case.name.replace("pml_", "enlarged_")
        s, alpha, rho, beta = media(case)
        _, alpha_ref, _, _ = media(ref)
        if np.abs(alpha - alpha_ref).max() > 1e-12:
            raise RuntimeError(f"{case}: drop differs from the reference")
        t, p = probes(case, "p")
        _, p_ref = probes(ref, "p")
        _, u = probes(case, "U")
        _, u_ref = probes(ref, "U")
        u = u.reshape(len(t), -1, 3)[:, :, :2]
        u_ref = u_ref.reshape(len(t), -1, 3)[:, :, :2]
        obs = t <= T_OBS * (1 + 1e-9)
        p, p_ref, u, u_ref = p[obs], p_ref[obs], u[obs], u_ref[obs]
        groups = {"all": np.ones_like(alpha, bool), "gas": alpha == 0,
                  "liquid": alpha == 1, "mixed": (alpha > 0) & (alpha < 1)}
        e = {g: e2(p, p_ref, beta, m) for g, m in groups.items()}
        ev = e2(u, u_ref, rho, groups["all"])
        row = [f"{s['h'] * 1e3:g}", int(groups["liquid"].sum()),
               int(groups["mixed"].sum())]
        for key, value in [("all", e["all"]), ("gas", e["gas"]),
                           ("liquid", e["liquid"]), ("mixed", e["mixed"]),
                           ("U", ev)]:
            q = (f"{math.log(prev[key] / value) / math.log(2):.2f}"
                 if prev else "")
            row += [f"{value:.3e}", q]
        prev = {"all": e["all"], "gas": e["gas"], "liquid": e["liquid"],
                "mixed": e["mixed"], "U": ev}
        rows.append(row)
    header = ["h_mm", "liquidCells", "mixedCells", "E2p", "q", "E2p_gas", "q",
              "E2p_liquid", "q", "E2p_mixed", "q", "E2U", "q"]
    write(root / "studyE_accuracy.csv", header, rows)


def long_time(root):
    rows = []
    for case in sorted((root / "long").iterdir()):
        path = case / "postProcessing/pmlVerification/summary.json"
        if not path.is_file():
            continue
        s = json.loads(path.read_text())
        e, m, q = s["energy"], s.get("maxima", {}), s["lastQuarter"]
        rows.append([case.name, f"{s['setup']['courant']:g}",
                     f"{e['atTobs']['Eac']:.3e}" if e.get("atTobs") else "",
                     f"{e['final']['Eac']:.3e}",
                     f"{m['p']['final']:.3e}", f"{m['psi']['final']:.3e}",
                     f"{m['w']['final']:.3e}", q["Eac"]["class"],
                     q["max_p"]["class"], q["max_psi"]["class"],
                     q["max_w"]["class"], f"{q['Eac']['rate']:.3e}",
                     f"{s['lastSampleTime']:.4g}"])
    write(root / "studyE_long.csv",
          ["case", "Co", "Eac_Tobs", "Eac_final", "maxP", "maxPsi", "maxW",
           "Eac", "maxP", "maxPsi", "maxW", "EacRate", "tLast"], rows)


def write(path, header, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"\n{path}")
    widths = [max(len(str(v)) for v in col) for col in zip(header, *rows)]
    for row in [header, *rows]:
        print("  ".join(str(v).rjust(w) for v, w in zip(row, widths)))


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "studies/E")
    if (root / "accuracy").is_dir():
        accuracy(root)
    if (root / "long").is_dir():
        long_time(root)


if __name__ == "__main__":
    main()
