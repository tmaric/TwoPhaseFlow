#!/usr/bin/env python3
"""Apply a boundary-preserving deformation to an ASCII OpenFOAM mesh."""

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
    parser.add_argument("--orthogonal-boundary-layers", type=int, default=0)
    parser.add_argument("--boundary-transition-fraction", type=float, default=0.0)
    args = parser.parse_args()

    if args.orthogonal_boundary_layers < 1:
        parser.error("--orthogonal-boundary-layers must be at least one")
    if args.boundary_transition_fraction < 0:
        parser.error("--boundary-transition-fraction must be nonnegative")

    path = pathlib.Path("constant/polyMesh/points")
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(
        r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)"
    )

    coordinates = [
        tuple(float(match.group(i)) for i in range(1, 4))
        for match in pattern.finditer(text)
    ]
    x_coordinates = sorted({coordinate[0] for coordinate in coordinates})
    y_coordinates = sorted({coordinate[1] for coordinate in coordinates})
    x_indices = {coordinate: index for index, coordinate in enumerate(x_coordinates)}
    y_indices = {coordinate: index for index, coordinate in enumerate(y_coordinates)}

    max_orthogonal_layers = (min(len(x_coordinates), len(y_coordinates)) - 1)//2
    if args.orthogonal_boundary_layers > max_orthogonal_layers:
        parser.error(
            "--orthogonal-boundary-layers exceeds half the smaller mesh dimension"
        )

    scale = args.amplitude * min(args.lx, args.ly)
    transition_width = (
        args.boundary_transition_fraction * min(args.lx, args.ly)
    )

    def axis_weight(coordinate: float, coordinates: list[float], index: int) -> float:
        upper_index = len(coordinates) - 1 - index
        distance_index = min(index, upper_index)
        if distance_index <= args.orthogonal_boundary_layers:
            return 0.0
        if transition_width == 0.0:
            return 1.0

        if index <= upper_index:
            boundary_distance = coordinate - coordinates[0]
            orthogonal_width = (
                coordinates[args.orthogonal_boundary_layers] - coordinates[0]
            )
        else:
            boundary_distance = coordinates[-1] - coordinate
            orthogonal_width = (
                coordinates[-1]
                - coordinates[-1 - args.orthogonal_boundary_layers]
            )

        normalized_distance = min(
            max((boundary_distance - orthogonal_width)/transition_width, 0.0),
            1.0,
        )
        return normalized_distance*normalized_distance*(3.0 - 2.0*normalized_distance)

    def deformation_weight(x: float, y: float) -> float:
        return min(
            axis_weight(x, x_coordinates, x_indices[x]),
            axis_weight(y, y_coordinates, y_indices[y]),
        )

    def deform(match: re.Match[str]) -> str:
        x, y, z = (float(match.group(i)) for i in range(1, 4))
        weight = deformation_weight(x, y)
        if weight == 0.0:
            return f"({x:.16g} {y:.16g} {z:.16g})"
        dx = weight*scale * math.sin(math.pi*x/args.lx) * math.sin(2*math.pi*y/args.ly)
        dy = weight*0.5*scale * math.sin(2*math.pi*x/args.lx) * math.sin(math.pi*y/args.ly)
        return f"({x + dx:.16g} {y + dy:.16g} {z:.16g})"

    new_text, count = pattern.subn(deform, text)
    if count == 0:
        raise RuntimeError(f"No mesh points found in {path}")
    path.write_text(new_text, encoding="utf-8")
    print(
        f"Deformed {count} mesh points with amplitude factor {args.amplitude}; "
        f"orthogonal boundary layers: {args.orthogonal_boundary_layers}; "
        f"transition fraction: {args.boundary_transition_fraction}"
    )


if __name__ == "__main__":
    main()
