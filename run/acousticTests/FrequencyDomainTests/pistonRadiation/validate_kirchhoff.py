#!/usr/bin/env python3
"""Validate the Kirchhoff quadrature with an exact outgoing spherical field."""

import math

import numpy as np

from postprocess_compare import (
    complete_sound_hard_meridian,
    kirchhoff_exterior_pressure,
    revolve_axisymmetric_source,
)


def main() -> int:
    frequency = 10000.0
    sound_speed = 343.0
    k = 2.0 * math.pi * frequency / sound_speed
    r_src = 0.2
    r_obs = 4.0

    theta = np.linspace(0.0, 0.5 * math.pi, 361)
    xyz = np.column_stack(
        [
            r_src * np.sin(theta),
            r_src * np.cos(theta),
            np.zeros_like(theta),
        ]
    )
    dtheta = np.gradient(theta)
    area = 2.0 * math.pi * r_src**2 * np.sin(theta) * dtheta
    pressure = np.full(theta.shape, np.exp(1j * k * r_src) / r_src)
    dpdn = (1j * k - 1.0 / r_src) * pressure

    valid = area > 0.0
    xyz, area, pressure, dpdn = complete_sound_hard_meridian(
        xyz[valid], area[valid], pressure[valid], dpdn[valid]
    )
    src_xyz, src_n, src_area, src_p, src_dpdn = revolve_axisymmetric_source(
        xyz, area, pressure, dpdn, 180
    )

    angles = np.deg2rad(np.array([0.0, 20.0, 45.0, 70.0, 90.0]))
    obs = np.column_stack(
        [
            r_obs * np.sin(angles),
            r_obs * np.cos(angles),
            np.zeros_like(angles),
        ]
    )
    reconstructed = kirchhoff_exterior_pressure(
        obs, src_xyz, src_n, src_area, src_p, src_dpdn, k
    )
    exact = np.full(angles.shape, np.exp(1j * k * r_obs) / r_obs)
    rel_l2 = float(np.linalg.norm(reconstructed - exact) / np.linalg.norm(exact))
    rel_linf = float(np.max(np.abs(reconstructed - exact)) / np.max(np.abs(exact)))

    print(f"complex relL2   = {rel_l2:.6e}")
    print(f"complex relLinf = {rel_linf:.6e}")
    if rel_linf > 2.0e-3:
        raise SystemExit("Kirchhoff reconstruction validation failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
