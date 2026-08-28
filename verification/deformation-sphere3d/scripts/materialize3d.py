"""Materialise one 3D spiralling-deformation case from the testsuite template.

Three files are touched: system/simulationParameter (resolution and
reconstruction scheme), system/controlDict (write interval and end time) and
system/fvSolution (the isoAdvector block carrying movingFrame3D).

The mesh is [0,1] x [0,1] x [0,2], so nz = 2 nx keeps the cells cubic.
"""

from __future__ import annotations

import os
import re
import shutil

_JUNK = ("volumeFractionError.dat", "results_deformation.csv", "error.dat")


def _iso_block(
    frame: str,
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
    rev = revolutions if frame == "frameA" else 0.0
    tra = translation_amplitude if frame == "frameA" else 0.0
    if frame not in ("none", "frameA"):
        raise ValueError(f"unknown frame: {frame!r}")

    return (
        "isoAdvector\n"
        "{\n"
        "    movingFrame3D\n"
        "    {\n"
        "        rotationCentre       %s;\n"
        "        swirlCentre          %s;\n"
        "        revolutions          %g;\n"
        "        endTime              %g;\n"
        "        period               %g;\n"
        "        baseAmplitude        1;\n"
        "        translationAmplitude %g;\n"
        "    }\n"
        "}"
    ) % (rotation_centre, swirl_centre, rev, end_time, period, tra)


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
    s = re.sub(r"^nz\s.*$", f"nz {2 * N};", s, flags=re.M)
    s = re.sub(r"^RECONSCHEME\s.*$", f"RECONSCHEME {scheme};", s, flags=re.M)
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
