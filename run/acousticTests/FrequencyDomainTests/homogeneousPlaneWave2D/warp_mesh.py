#!/usr/bin/env python3
"""Apply a smooth boundary-preserving deformation to an ASCII OpenFOAM mesh."""

from __future__ import annotations

import argparse
import math
import pathlib
import re


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--amplitude", type=float, required=True)
    parser.add_argument("--lx", type=float, required=True)
    parser.add_argument("--ly", type=float, required=True)
    args = parser.parse_args()

    path = pathlib.Path("constant/polyMesh/points")
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(
        r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)"
    )

    scale = args.amplitude * min(args.lx, args.ly)

    def deform(match: re.Match[str]) -> str:
        x, y, z = (float(match.group(i)) for i in range(1, 4))
        dx = scale * math.sin(math.pi*x/args.lx) * math.sin(2*math.pi*y/args.ly)
        dy = 0.5*scale * math.sin(2*math.pi*x/args.lx) * math.sin(math.pi*y/args.ly)
        return f"({x + dx:.16g} {y + dy:.16g} {z:.16g})"

    new_text, count = pattern.subn(deform, text)
    if count == 0:
        raise RuntimeError(f"No mesh points found in {path}")
    path.write_text(new_text, encoding="utf-8")
    print(f"Warped {count} mesh points with amplitude factor {args.amplitude}")


if __name__ == "__main__":
    main()

