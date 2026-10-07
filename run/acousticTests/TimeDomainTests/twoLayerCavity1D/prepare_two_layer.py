#!/usr/bin/env python3
"""Generate the run-dependent inputs of the two-layer cavity from the
environment set by prepareCase."""
import json
import math
import os
from pathlib import Path

from twoLayerModes import TwoLayerMode, modes_near, select_mode

CASE = Path(__file__).resolve().parent


def env(name, kind=float):
    return kind(os.environ[name])


def fill(template, target, values):
    text = (CASE / template).read_text()
    for key, value in values.items():
        text = text.replace(f"__{key}__", str(value))
    (CASE / target).write_text(text)


def main():
    rho_g, c_g, rho_l, c_l = (env(n) for n in ("RHOG", "CG", "RHOL", "CL"))
    gas_units, liquid_units = env("GAS_UNITS", int), env("LIQUID_UNITS", int)
    per_unit = env("CELLS_PER_UNIT", int)
    unit = env("L_GAS") / gas_units
    h = unit / per_unit
    n_gas, n_liquid = gas_units * per_unit, liquid_units * per_unit
    nx = n_gas + n_liquid
    length = nx * h
    offset = env("INTERFACE_OFFSET")
    if not 0 <= offset < 1:
        raise SystemExit("INTERFACE_OFFSET must lie in [0, 1)")
    x_i = n_gas * h + offset * h

    roots = modes_near(env("F_TARGET"), rho_g, c_g, rho_l, c_l, x_i, length)
    omega = select_mode(roots, env("F_TARGET"), os.environ["MODE"])
    mode = TwoLayerMode(rho_g, c_g, rho_l, c_l, x_i, length, omega,
                        env("PA"))
    freq = omega / (2 * math.pi)

    steps_per_period = env("STEPS_PER_PERIOD", int)
    samples = env("PROBE_SAMPLES_PER_PERIOD", int)
    stride = 1 if steps_per_period <= samples else steps_per_period // samples
    if stride * samples != steps_per_period and stride != 1:
        raise SystemExit("STEPS_PER_PERIOD must be a multiple of "
                         "PROBE_SAMPLES_PER_PERIOD")
    n_steps = steps_per_period * env("N_PERIODS", int)
    delta_t = 1.0 / (freq * steps_per_period)
    end_time = n_steps * delta_t
    nprocs = env("NPROCS", int)

    fill("constant/transportProperties.in", "constant/transportProperties", {
        "RHOG": rho_g, "CG": c_g, "RHOL": rho_l, "CL": c_l,
        "F": repr(freq), "TIME_METHOD": os.environ["TIME_METHOD"],
        "THETA": os.environ["THETA"],
        "BLOCK_SOLVER": os.environ["BLOCK_SOLVER"],
    })
    fill("system/blockMeshDict.in", "system/blockMeshDict", {
        "L": repr(length), "HALF_H": repr(0.5 * h), "NX": nx,
    })
    fill("system/controlDict.in", "system/controlDict", {
        "SOLVER": os.environ["SOLVER"], "DELTA_T": repr(delta_t),
        "END_TIME": repr(end_time), "PROBE_STRIDE": stride,
    })
    fill("system/setExprFieldsDict.in", "system/setExprFieldsDict", {
        "XI": repr(x_i), "H": repr(h), "L": repr(length),
        "A": repr(mode.a), "B": repr(mode.b),
        "KG": repr(mode.k_g), "KL": repr(mode.k_l),
    })
    fill("system/decomposeParDict.in", "system/decomposeParDict",
         {"NPROCS": nprocs})

    centres = [(i + 0.5) * h for i in range(nx)]
    (CASE / "system/probeLocations").write_text(
        "probeLocations\n(\n"
        + "".join(f"    ({x!r} 0 0)\n" for x in centres)
        + ");\n"
    )

    if nprocs == 2:
        # Rank 0 holds the pure gas cells; with an offset the interface
        # cell opens rank 1, so it is adjacent to the processor face.
        owners = [0 if i < n_gas else 1 for i in range(nx)]
        (CASE / "constant/cellDecomposition").write_text(
            "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
            "    class labelList;\n    object cellDecomposition;\n}\n\n"
            f"{nx}\n(\n" + "\n".join(map(str, owners)) + "\n)\n"
        )
    elif nprocs != 1:
        raise SystemExit("NPROCS must be 1 or 2")

    setup = {
        "rhoG": rho_g, "cG": c_g, "rhoL": rho_l, "cL": c_l,
        "length": length, "xInterface": x_i, "interfaceOffset": offset,
        "nGas": n_gas, "nLiquid": n_liquid, "nx": nx, "h": h,
        "crossSection": h * h, "cellsPerUnit": per_unit,
        "gasCellsPerWavelength": c_g / freq / h,
        "mode": os.environ["MODE"], "omega": omega, "f": freq,
        "modesNearTargetHz": [w / (2 * math.pi) for w in roots],
        "A": mode.a, "B": mode.b, "Pa": env("PA"),
        "deltaT": delta_t, "nSteps": n_steps,
        "stepsPerPeriod": steps_per_period, "nPeriods": env("N_PERIODS", int),
        "probeStride": stride, "method": os.environ["TIME_METHOD"],
        "theta": env("THETA"), "blockSolver": os.environ["BLOCK_SOLVER"],
        "nProcs": nprocs,
    }
    (CASE / "caseSetup.json").write_text(json.dumps(setup, indent=2) + "\n")
    print(f"mode {setup['mode']}: f = {freq:.6f} Hz, nx = {nx}, "
          f"x_I = {x_i:.9g} m, dt = {delta_t:.6g} s")


if __name__ == "__main__":
    main()
