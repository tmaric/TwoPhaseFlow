"""Materialise one reversed-vortex case from the testsuite template.

Only three files are touched: system/simulationParameter (resolution and
reconstruction scheme), system/controlDict (write interval and end time) and
system/fvSolution (the isoAdvector block, which carries the movingFrame
sub-dictionary when a frame is requested).
"""

from __future__ import annotations

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
                        u_max: float) -> None:
    """Write maxAlphaCo, and an initial step that respects it.

    The template carries a fixed deltaT, so the FIRST step runs at an interface
    Courant number of roughly deltaT * N whatever the mesh -- at N = 2048 that
    was 2.13 against a limit of 0.5. Above Courant 1 the swept volume can reach
    past the neighbouring cell and the geometric flux construction is outside
    its design range, and the error that injects does not heal. Scaling the
    initial step with both the mesh and the limit keeps the first step well
    inside it whatever the field's magnitude, and the controller takes over.
    """
    ctrl = os.path.join(case, "system", "controlDict")
    with open(ctrl) as fh:
        s = fh.read()
    s = re.sub(r"^maxAlphaCo\s.*$", f"maxAlphaCo      {max_alpha_co:g};", s, flags=re.M)
    # no constant cap: it would assume |u| of order one, and the LeVeque
    # field peaks near 2. 0.25*co/N keeps the first step at or below
    # 0.1 for |u| <= 2, and the controller reaches the target in a few steps.
    dt = 0.25 * max_alpha_co / (float(N) * u_max)
    # The ceiling matters as much as the initial step: it is what the
    # controller grows to while the base amplitude passes through zero,
    # and a mesh-independent ceiling turns into an unbounded Courant
    # number under refinement.
    max_dt = max_alpha_co / (float(N) * u_max)
    s = re.sub(r"^maxDeltaT\s.*$", f"maxDeltaT       {max_dt:g};", s, flags=re.M)
    s = re.sub(r"^deltaT\s.*$", f"deltaT          {dt:g};", s, flags=re.M)
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
    u_max: float = 1.05,
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
    _write_time_control(case, N, max_alpha_co, u_max)
