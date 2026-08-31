"""Materialise one reversed-vortex case from the testsuite template.

Only three files are touched: system/simulationParameter (resolution and
reconstruction scheme), system/controlDict (write interval and end time) and
system/fvSolution (the isoAdvector block, which carries the movingFrame
sub-dictionary when a frame is requested).
"""

from __future__ import annotations

import math
import os
import re
import shutil

# Files that are outputs of a previous run and must never be inherited.
_JUNK = (
    "volumeFractionError.dat",
    "volumeFractionErrorPar.dat",
    "error.dat",
    "results_shearedDisc.csv",
)


def _iso_block(
    frame: str,
    period: float,
    end_time: float,
    rotation_centre: str,
    revolutions: float,
    translation_amplitude: float,
) -> str:
    """The isoAdvector dictionary, with or without the moving frame."""
    if frame == "none":
        return "isoAdvector\n{\n    period      %g;\n}" % period
    if frame == "frameId":
        # The identity frame: the moving-frame code path with the frame itself
        # switched off, so the result must reproduce the unframed run exactly.
        revolutions = 0.0
        translation_amplitude = 0.0
    elif frame != "frameA":
        raise ValueError(f"unknown frame: {frame!r}")
    return (
        "isoAdvector\n"
        "{\n"
        "    period      %g;\n"
        "    movingFrame\n"
        "    {\n"
        "        rotationCentre       %s;\n"
        "        revolutions          %g;\n"
        "        endTime              %g;\n"
        "        period               %g;\n"
        "        baseAmplitude        1;\n"
        "        translationAmplitude %g;\n"
        "    }\n"
        "}"
    ) % (
        period,
        rotation_centre,
        revolutions,
        end_time,
        period,
        translation_amplitude,
    )


def _write_time_control(case: str, N: int, max_alpha_co: float,
                        u_max: float, end_time: float) -> None:
    """Fix the time step from the known maximum speed.

    dt = Co*h/max|w|, with max|w| the exact maximum of the laboratory field
    over the domain and the period. The temporal factor has modulus at most
    one, so this is known before the run and needs no controller: the Courant
    number is Co in every cell at every step.

    The step count is rounded up to a multiple of four so that t = T/4, T/2,
    3T/4 and T all land exactly on a step and the write-time adjustment is a
    no-op -- OpenFOAM's adjustableRunTime otherwise shifts dt to reach a write
    time, which makes the write interval part of the protocol.
    """
    ctrl = os.path.join(case, "system", "controlDict")
    with open(ctrl) as fh:
        s = fh.read()
    # u_max is max||w||_1, because OpenFOAM's Courant number is
    # 0.5*sum_f|phi_f|*dt/V = (|wx|+|wy|+|wz|)*dt/h on a Cartesian cell
    dt_target = max_alpha_co / (float(N) * u_max)
    steps = 4 * int(math.ceil(end_time / (4.0 * dt_target)))
    dt = end_time / steps
    s = re.sub(r"^adjustTimeStep\s.*$", "adjustTimeStep  no;", s, flags=re.M)
    s = re.sub(r"^deltaT\s.*$", f"deltaT          {dt:.10g};", s, flags=re.M)
    s = re.sub(r"^maxAlphaCo\s.*$", f"maxAlphaCo      {max_alpha_co:g};", s, flags=re.M)
    with open(ctrl, "w") as fh:
        fh.write(s)


def _write_decomposition(case: str, np: int) -> None:
    """Size system/decomposeParDict for this run.

    The templates ship a hard-coded `numberOfSubdomains 4`, and the two
    three-dimensional ones use `method simple`, which additionally needs an
    explicit `n (x y z)` and therefore does not survive an arbitrary rank count.
    scotch needs only the subdomain count.
    """
    p = os.path.join(case, "system", "decomposeParDict")
    if not os.path.exists(p):
        return
    s = open(p).read()
    s = re.sub(r"^numberOfSubdomains\s.*$", f"numberOfSubdomains {np};", s, flags=re.M)
    s = re.sub(r"^method\s.*$", "method          scotch;", s, flags=re.M)
    open(p, "w").write(s)


def build_case(
    template: str,
    case: str,
    scheme: str,
    frame: str,
    N: int,
    end_time: float,
    period: float,
    write_interval: float,
    rotation_centre: str,
    revolutions: float,
    translation_amplitude: float,
    np: int = 1,
    max_alpha_co: float = 0.2,
    u_max: float = 1.607679,
) -> None:
    if os.path.exists(case):
        shutil.rmtree(case)
    shutil.copytree(template, case)

    for name in _JUNK:
        p = os.path.join(case, name)
        if os.path.exists(p):
            os.remove(p)

    # resolution + reconstruction scheme
    p = os.path.join(case, "system", "simulationParameter")
    s = open(p).read()
    s = re.sub(r"^nx\s.*$", f"nx\t{N};", s, flags=re.M)
    s = re.sub(r"^nz\s.*$", f"nz\t{N};", s, flags=re.M)
    s = re.sub(r"^RECONSCHEME\s.*$", f"RECONSCHEME {scheme};", s, flags=re.M)
    open(p, "w").write(s)

    # write interval (top level and inside the function object) + end time
    p = os.path.join(case, "system", "controlDict")
    s = open(p).read()
    s = re.sub(r"^(\s*)writeInterval\s.*$", rf"\g<1>writeInterval   {write_interval:g};", s, flags=re.M)
    s = re.sub(r"^endTime\s.*$", f"endTime         {end_time:g};", s, flags=re.M)
    # keep the fields human-readable so foamlib can parse them without gunzip
    s = re.sub(r"^writeCompression\s.*$", "writeCompression off;", s, flags=re.M)
    open(p, "w").write(s)

    # isoAdvector block
    p = os.path.join(case, "system", "fvSolution")
    s = open(p).read()
    block = _iso_block(
        frame, period, end_time, rotation_centre, revolutions, translation_amplitude
    )
    s, n = re.subn(r"isoAdvector\s*\{[^}]*\}", lambda _m: block, s, count=1)
    if n != 1:
        raise RuntimeError(f"could not find the isoAdvector block in {p}")
    open(p, "w").write(s)

    _write_decomposition(case, np)
    _write_time_control(case, N, max_alpha_co, u_max, end_time)
