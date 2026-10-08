#!/usr/bin/env python3
"""Generate the run-dependent inputs of the piston-pulse case from the
environment set by prepareCase. Run after 0.orig has been restored to 0."""
import json
import os
from pathlib import Path

from pulseReference import LayeredPulseSolution, PistonPulse

CASE = Path(__file__).resolve().parent


def env(name, kind=float):
    return kind(os.environ[name])


def fill(template, target, values):
    text = (CASE / template).read_text()
    for key, value in values.items():
        text = text.replace(f"__{key}__", str(value))
    (CASE / target).write_text(text)


def main():
    config = os.environ["CONFIG"]
    gas = (env("RHOG"), env("CG"))
    liquid = (env("RHOL"), env("CL"))
    if config == "gasFirst":
        (rho1, c1), (rho2, c2) = gas, liquid
    elif config == "liquidFirst":
        (rho1, c1), (rho2, c2) = liquid, gas
    else:
        raise SystemExit("CONFIG must be gasFirst or liquidFirst")

    length, x_i = env("LENGTH"), env("X_INTERFACE")
    nx = env("NX", int)
    h = length / nx
    if abs(x_i / h - round(x_i / h)) > 1e-9:
        raise SystemExit(f"NX = {nx} does not place X_INTERFACE on a face")
    pml_start = length - env("PML_THICKNESS")
    if pml_start <= x_i:
        raise SystemExit("the PML must lie in medium 2")

    f0 = env("F0")
    pulse = PistonPulse(f0, env("N_CYCLES"), env("U0"))
    solution = LayeredPulseSolution(pulse, rho1, c1, rho2, c2, x_i)

    steps_per_period = env("STEPS_PER_PERIOD", int)
    samples = env("PROBE_SAMPLES_PER_PERIOD", int)
    if steps_per_period % samples:
        raise SystemExit("STEPS_PER_PERIOD must be a multiple of "
                         "PROBE_SAMPLES_PER_PERIOD")
    stride = steps_per_period // samples
    delta_t = 1.0 / (f0 * steps_per_period)
    n_steps = round(env("T_OBS") * f0 * steps_per_period)
    end_time = n_steps * delta_t

    fill("constant/transportProperties.in", "constant/transportProperties", {
        "RHOG": gas[0], "CG": gas[1], "RHOL": liquid[0], "CL": liquid[1],
        "F0": f0, "TIME_METHOD": os.environ["TIME_METHOD"],
        "THETA": os.environ["THETA"], "BLOCK_SOLVER": os.environ["BLOCK_SOLVER"],
        "PML_THICKNESS": os.environ["PML_THICKNESS"],
        "SIGMA_MAX": os.environ["SIGMA_MAX"],
        "PML_ORDER": os.environ["PML_ORDER"], "LENGTH": length,
    })
    fill("system/blockMeshDict.in", "system/blockMeshDict", {
        "L": repr(length), "HALF_H": repr(0.5 * h), "NX": nx,
    })
    fill("system/controlDict.in", "system/controlDict", {
        "SOLVER": os.environ["SOLVER"], "DELTA_T": repr(delta_t),
        "END_TIME": repr(end_time), "PROBE_STRIDE": stride,
    })
    liquid_first = config == "liquidFirst"
    fill("system/setExprFieldsDict.in", "system/setExprFieldsDict", {
        "LIQUID_SIDE": "left" if liquid_first else "right",
        "ALPHA_EXPRESSION":
            f"(pos().x() {'<' if liquid_first else '>'} {x_i!r}) ? 1 : 0",
    })
    fill("0.orig/p.in", "0/p", {
        "DRIVE_EXPRESSION": pulse.acceleration_expression(rho1),
    })
    (CASE / "0" / "p.in").unlink(missing_ok=True)

    (CASE / "system/probeLocations").write_text(
        "probeLocations\n(\n"
        + "".join(f"    ({(i + 0.5) * h!r} 0 0)\n" for i in range(nx))
        + ");\n"
    )

    setup = {
        "config": config, "rho1": rho1, "c1": c1, "rho2": rho2, "c2": c2,
        "reflection": solution.r, "length": length, "xInterface": x_i,
        "pmlStart": pml_start, "pmlThickness": env("PML_THICKNESS"),
        "sigmaMax": env("SIGMA_MAX"), "pmlOrder": env("PML_ORDER"),
        "nx": nx, "h": h, "gasCellsPerWavelength": gas[1] / f0 / h,
        "f0": f0, "nCycles": env("N_CYCLES"), "u0": env("U0"),
        "deltaT": delta_t, "nSteps": n_steps, "endTime": end_time,
        "stepsPerPeriod": steps_per_period, "probeStride": stride,
        "method": os.environ["TIME_METHOD"], "theta": env("THETA"),
        "blockSolver": os.environ["BLOCK_SOLVER"],
    }
    (CASE / "caseSetup.json").write_text(json.dumps(setup, indent=2) + "\n")
    print(f"{config}: nx = {nx}, N_lambda_g = {setup['gasCellsPerWavelength']:.1f}, "
          f"dt = {delta_t:.4g} s, {n_steps} steps, R = {solution.r:+.6f}")


if __name__ == "__main__":
    main()
