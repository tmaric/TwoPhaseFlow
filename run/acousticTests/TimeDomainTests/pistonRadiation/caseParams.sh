#!/bin/sh

# Single source of truth for the time-domain pistonRadiation case.
# Environment variables may override the more expensive run controls.

SOLVER=acousticWaveFoam

# Drive and medium
DRIVE_F=${DRIVE_F:-10000}
PISTON_U=${PISTON_U:-1}
RHOL=998.3
RHOG=1.2
CG=343
CL=1480

# Geometry and radial PML
PISTON_RADIUS=0.1
PML_RMIN=0.2
PML_L=0.08

# Mesh
MESH_CELLS_PER_LAMBDA=${MESH_CELLS_PER_LAMBDA:-60}
MESH_GR=1.0
WEDGE_ANGLE_DEG=2

# Time integration and output
N_PERIODS=${N_PERIODS:-25}
PROBE_SAMPLES_PER_PERIOD=${PROBE_SAMPLES_PER_PERIOD:-50}
