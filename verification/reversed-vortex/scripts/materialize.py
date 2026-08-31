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


def _write_initial_deltat(case: str, N: int) -> None:
    """Keep the first step inside maxAlphaCo.

    The template's fixed deltaT = 1e-3 gives a first-step interface Courant
    number of about 1e-3 * N, so from N = 512 upwards the very first step
    violates the maxAlphaCo = 0.5 the study runs under, and above Courant 1
    the geometric swept-volume assumption fails outright. The controller
    recovers on step two, but the volume and boundedness error it injects
    does not. Levels at or below 256 are left untouched, so their agreement
    with the published values is preserved bit for bit.
    """
    ctrl = os.path.join(case, "system", "controlDict")
    with open(ctrl) as fh:
        s = fh.read()
    dt = min(1.0e-3, 0.4 / float(N))
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
    _write_initial_deltat(case, N)
