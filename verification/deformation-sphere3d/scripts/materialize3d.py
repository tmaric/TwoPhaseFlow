"""Materialise one 3D advection case from a testsuite template.

Three files are touched: system/simulationParameter (resolution and
reconstruction scheme), system/controlDict (write interval and end time) and
system/fvSolution (the isoAdvector block carrying movingFrame3D).

`nz_factor * nx` keeps the cells cubic: the spiralling case lives on
[0,1] x [0,1] x [0,2] (factor 2), the LeVeque case on [0,1]^3 (factor 1).
Note that the committed simulationParameter of test-Leveque carries a stale
`nz 64` against a unit-cube domain; the testsuite overrides it with nz = res,
and so does this workflow.
"""

from __future__ import annotations

import math
import os
import re
import shutil

_JUNK = (
    "volumeFractionError.dat",
    "results_deformation.csv",
    "results_leveque.csv",
    "error.dat",
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


def _iso_block(
    frame: str,
    base_field: str,
    period: float,
    end_time: float,
    rotation_centre: str,
    swirl_centre: str,
    revolutions: float,
    translation_amplitude: float,
) -> str:
    """The isoAdvector dictionary.

    Both columns go through movingFrameFlow3D, so the only difference between
    them is the frame itself: `none` sets zero revolutions and zero
    translation, which leaves the plain forward-backward base field.
    """
    if frame not in ("none", "frameA"):
        raise ValueError(f"unknown frame: {frame!r}")
    rev = revolutions if frame == "frameA" else 0.0
    tra = translation_amplitude if frame == "frameA" else 0.0

    lines = [
        "isoAdvector",
        "{",
        "    movingFrame3D",
        "    {",
        f"        baseField            {base_field};",
        f"        rotationCentre       {rotation_centre};",
        f"        swirlCentre          {swirl_centre};",
        f"        revolutions          {rev:g};",
        f"        endTime              {end_time:g};",
        f"        period               {period:g};",
        "        baseAmplitude        1;",
        f"        translationAmplitude {tra:g};",
        "    }",
        "}",
    ]
    return "\n".join(lines)


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
    swirl_centre: str,
    revolutions: float,
    translation_amplitude: float,
    base_field: str = "spiralling",
    nz_factor: int = 2,
    np: int = 1,
    max_alpha_co: float = 0.2,
    u_max: float = 2.344226,
    domain_height: float = 2.0,
) -> None:
    if os.path.exists(case):
        shutil.rmtree(case)
    shutil.copytree(template, case)

    for name in _JUNK:
        p = os.path.join(case, name)
        if os.path.exists(p):
            os.remove(p)

    p = os.path.join(case, "system", "simulationParameter")
    s = open(p).read()
    s = re.sub(r"^nx\s.*$", f"nx {N};", s, flags=re.M)
    s = re.sub(r"^ny\s.*$", f"ny {N};", s, flags=re.M)
    s = re.sub(r"^nz\s.*$", f"nz {nz_factor * N};", s, flags=re.M)
    s = re.sub(r"^RECONSCHEME\s.*$", f"RECONSCHEME {scheme};", s, flags=re.M)
    open(p, "w").write(s)

    # domain height: the spiralling case ships a box of height 2, but the
    # material only reaches z = 0.87, so the upper half is dead mesh
    p = os.path.join(case, "system", "blockMeshDict")
    s = open(p).read()
    s = re.sub(r"^z2\s.*$", f"z2 {domain_height:g};", s, flags=re.M)
    open(p, "w").write(s)

    p = os.path.join(case, "system", "controlDict")
    s = open(p).read()
    s = re.sub(
        r"^(\s*)writeInterval\s.*$",
        rf"\g<1>writeInterval   {write_interval:g};",
        s,
        flags=re.M,
    )
    s = re.sub(r"^endTime\s.*$", f"endTime         {end_time:g};", s, flags=re.M)
    # readable fields, so the volume fraction can be parsed without gunzip
    if re.search(r"^writeCompression\s", s, flags=re.M):
        s = re.sub(r"^writeCompression\s.*$", "writeCompression off;", s, flags=re.M)
    else:
        s = s.replace("writeFormat", "writeCompression off;\n\nwriteFormat", 1)
    open(p, "w").write(s)

    p = os.path.join(case, "system", "fvSolution")
    s = open(p).read()
    block = _iso_block(
        frame,
        base_field,
        period,
        end_time,
        rotation_centre,
        swirl_centre,
        revolutions,
        translation_amplitude,
    )
    s, n = re.subn(r"isoAdvector\s*\{[^}]*\}", lambda _m: block, s, count=1)
    if n != 1:
        raise RuntimeError(f"could not find the isoAdvector block in {p}")
    open(p, "w").write(s)

    _write_decomposition(case, np)
    _write_time_control(case, N, max_alpha_co, u_max, end_time)
