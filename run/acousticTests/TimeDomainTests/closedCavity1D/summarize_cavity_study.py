#!/usr/bin/env python3
"""Tabulate the closed-cavity studies A, B and C written by runStudy."""
import csv
import json
import math
import sys
from pathlib import Path


def load(study_dir):
    rows = []
    for path in sorted(study_dir.glob("*/postProcessing/cavityVerification/"
                                      "summary.json")):
        data = json.loads(path.read_text())
        data["name"] = path.parents[2].name
        rows.append(data)
    return rows


def orders(values, ratio=2.0):
    out = [None]
    for coarse, fine in zip(values, values[1:]):
        ok = coarse and fine and coarse > 0 and fine > 0
        out.append(math.log(coarse / fine) / math.log(ratio) if ok else None)
    return out


def write(path, header, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"\n{path}")
    widths = [max(len(str(v)) for v in col) for col in zip(header, *rows)]
    for row in [header, *rows]:
        print("  ".join(str(v).rjust(w) for v, w in zip(row, widths)))


def fmt(value, digits=4):
    return "" if value is None else f"{value:.{digits}e}"


def study_a(out, rows):
    table = []
    for r in rows:
        rates = r["rates"]
        table.append([
            r["name"], r["setup"]["method"], r["setup"]["theta"],
            r["setup"]["blockSolver"],
            fmt(r["implementation"]["maxModalErrorVsDiscrete"], 3),
            fmt(r["E2"]["vsDiscrete"], 3),
            fmt(rates["gammaFit"]), fmt(rates["gammaDiscrete"]),
            fmt(rates["omegaFit"] / rates["omegaDiscrete"] - 1.0, 3),
            fmt(r["energy"]["discreteWaveEnergyMaxRelativeDrift"], 3),
            fmt(r["energy"]["physicalEnergyLastPeriodRelativeSwing"], 3),
        ])
    write(out / "studyA.csv",
          ["case", "method", "theta", "solver", "maxModalErrDisc",
           "E2disc", "gammaFit", "gammaDisc", "omegaFit/omegaDisc-1",
           "EhDrift", "EacSwing"], table)


def study_b(out, rows):
    for method in ("blockCoupled", "fullyImplicit"):
        sel = sorted((r for r in rows if r["setup"]["method"] == method),
                     key=lambda r: r["setup"]["stepsPerPeriod"])
        if not sel:
            continue
        e2 = [r["E2"]["vsSemiDiscrete"] for r in sel]
        freq = [abs(r["rates"]["relativeFrequencyErrorVsSemiDiscrete"])
                for r in sel]
        damp = [r["rates"]["gammaFit"] / r["omega_h"] for r in sel]
        o_e2, o_freq, o_damp = orders(e2), orders(freq), orders(damp)
        table = [[
            r["setup"]["stepsPerPeriod"], f"{r['xi']:.4f}",
            f"{r['courant']:.3g}", fmt(e2[i]), fmt(o_e2[i], 2),
            fmt(freq[i]), fmt(o_freq[i], 2), fmt(damp[i]), fmt(o_damp[i], 2),
            fmt(r["endOfRun"]["amplitudeFactorFit"]),
        ] for i, r in enumerate(sel)]
        write(out / f"studyB_{method}.csv",
              ["N_T", "xi", "Courant", "E2semi", "order", "|freqErr|",
               "order", "gamma/omega_h", "order", "ampFactor"], table)


def study_c(out, rows):
    sel = sorted(rows, key=lambda r: r["setup"]["cellsPerWavelength"])
    spatial = [abs(r["rates"]["relativeSpatialFrequencyError"]) for r in sel]
    e2 = [r["E2"]["vsExact"] for r in sel]
    e2u = [r["E2"]["velocityVsExact"] for r in sel]
    o_sp, o_e2, o_e2u = orders(spatial), orders(e2), orders(e2u)
    table = [[
        r["setup"]["cellsPerWavelength"], r["setup"]["nx"],
        f"{r['setup']['k'] * r['setup']['h']:.4f}", fmt(spatial[i]),
        fmt(abs(r["rates"]["relativeSpatialFrequencyErrorPredicted"])),
        fmt(o_sp[i], 2), fmt(e2[i]), fmt(o_e2[i], 2), fmt(e2u[i]),
        fmt(o_e2u[i], 2),
        fmt(r["energy"]["physicalEnergyLastPeriodRelativeSwing"], 3),
    ] for i, r in enumerate(sel)]
    write(out / "studyC.csv",
          ["N_lambda", "nx", "kh", "|omega_h/omega-1|", "predicted",
           "order", "E2p", "order", "E2U", "order", "EacSwing"], table)


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "studies")
    for key, handler in (("A", study_a), ("B", study_b), ("C", study_c)):
        rows = load(out / key)
        if rows:
            handler(out, rows)


if __name__ == "__main__":
    main()
