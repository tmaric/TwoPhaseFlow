#!/usr/bin/env python3
"""Tabulate the multidirectional PML studies written by runStudy.

Study A: E2 against the exact free-space solution for the PML and the
enlarged-domain runs, and the PML truncation error, E2 of the PML run
against the enlarged-domain run on the same cells, time-integrated over
the samples up to T_obs.
Study B: physical-domain energy and whole-domain maxima, with their
last-quarter classification.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

T_OBS = 4e-4


def summary(case):
    return json.loads((case / "postProcessing/pmlVerification/"
                       "summary.json").read_text())


def pressure(case):
    rows = [line.split() for line in
            (case / "postProcessing/cellProbes/0/p").read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    data = np.array(rows, dtype=float)
    return data[:, 0], data[:, 1:]


def fmt(x, digits=3):
    return "" if x is None else f"{x:.{digits}e}"


def order(coarse, fine):
    if not coarse or not fine:
        return ""
    return f"{math.log(coarse / fine) / math.log(2.0):.2f}"


def write(path, header, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"\n{path}")
    widths = [max(len(str(v)) for v in col) for col in zip(header, *rows)]
    for row in [header, *rows]:
        print("  ".join(str(v).rjust(w) for v, w in zip(row, widths)))


def study_a(root):
    rows, previous = [], None
    for case in sorted((root / "A").glob("pml_h*"),
                       key=lambda c: -float(c.name[5:])):
        reference = case.parent / case.name.replace("pml_", "enlarged_")
        s, r = summary(case), summary(reference)
        t, p = pressure(case)
        t_ref, p_ref = pressure(reference)
        if not np.allclose(t, t_ref) or p.shape != p_ref.shape:
            raise RuntimeError(f"{case} and {reference} sample differently")
        obs = t <= T_OBS * (1 + 1e-9)
        # Homogeneous gas: the energy weight V beta is a constant factor.
        diff = p[obs] - p_ref[obs]
        truncation = math.sqrt(np.sum(diff**2) / np.sum(p_ref[obs] ** 2))
        peak = float(np.abs(diff).max() / np.abs(p_ref[obs]).max())
        e_pml = s["E2exact"]["pressure"]
        e_ref = r["E2exact"]["pressure"]
        rows.append([
            f"{s['setup']['h'] * 1e3:g}",
            f"{s['setup']['pulseWidth'] / s['setup']['h']:g}",
            fmt(e_ref), order(previous[0], e_ref) if previous else "",
            fmt(e_pml), order(previous[1], e_pml) if previous else "",
            fmt(r["E2exact"].get("velocity")), fmt(truncation), fmt(peak),
        ])
        previous = (e_ref, e_pml)
    write(root / "studyA.csv",
          ["h_mm", "s_over_h", "E2p_enlarged", "order", "E2p_pml", "order",
           "E2U_enlarged", "E2p_pml_vs_enlarged", "maxRel_pml_vs_enlarged"],
          rows)


def study_b(root):
    rows = []
    items = [summary(c) for c in (root / "B").iterdir()
             if (c / "postProcessing/pmlVerification/summary.json").is_file()]
    for s in sorted(items, key=lambda s: (
            s["setup"]["method"], s["setup"]["pmlThickness"],
            s["setup"]["courant"], s["setup"]["sigmaMax"])):
        st, e, q = s["setup"], s["energy"], s["lastQuarter"]
        m = s.get("maxima", {})
        rows.append([
            "monolithic" if st["method"] == "blockCoupled" else "segregated",
            f"{st['pmlThickness']:g}", f"{st['courant']:g}",
            f"{st['sigmaMax']:.3g}", f"{st['sigmaMax'] * st['deltaT']:.3g}",
            fmt(e["atTobs"]["Eac"]) if e.get("atTobs") else "",
            fmt(e["final"]["Eac"]),
            fmt(m.get("p", {}).get("final")),
            fmt(m.get("psi", {}).get("final")),
            fmt(m.get("w", {}).get("final")),
            q["Eac"]["class"], q.get("max_p", {}).get("class", ""),
            q.get("max_psi", {}).get("class", ""),
            q.get("max_w", {}).get("class", ""),
            fmt(s.get("blowupTime")), f"{s['lastSampleTime']:.4g}",
        ])
    write(root / "studyB.csv",
          ["update", "L", "Co", "sigmaMax", "sigma_dt", "Eac_Tobs",
           "Eac_final", "maxP_final", "maxPsi_final", "maxW_final",
           "Eac", "maxP", "maxPsi", "maxW", "blowup", "tLast"], rows)


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "studies")
    if (root / "A").is_dir():
        study_a(root)
    if (root / "B").is_dir():
        study_b(root)


if __name__ == "__main__":
    main()
