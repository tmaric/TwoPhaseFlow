#!/bin/sh

SOLVER=acousticHelmholtzFoam

# Nondimensional homogeneous medium. With c=1 m/s and f=1 Hz, lambda=1 m.
DRIVE_F=${DRIVE_F:-1}
RHO0=1
C0=1
P0=1
PHASE=0.5235987755982988

# Rectangular two-dimensional domain.
LX=1.125
LY=0.625
Z_HALF=0.005

# Study controls, overridable from the environment.
CELLS_PER_WAVELENGTH=${CELLS_PER_WAVELENGTH:-32}
# Supported mesh families:
#   orthogonal     - Cartesian mesh
#   warpedInterior - orthogonal boundary-cell layer with a non-orthogonal interior
MESH_FAMILY=${MESH_FAMILY:-warpedInterior}
# Supported boundary modes: dirichlet, alignedDirichlet, mixed.
BOUNDARY_MODE=${BOUNDARY_MODE:-dirichlet}
WARP_AMPLITUDE=${WARP_AMPLITUDE:-0.08}
ORTHOGONAL_BOUNDARY_LAYERS=${ORTHOGONAL_BOUNDARY_LAYERS:-1}
BOUNDARY_TRANSITION_FRACTION=${BOUNDARY_TRANSITION_FRACTION:-0.15}
N_NONORTH_CORRECTORS=${N_NONORTH_CORRECTORS:-2}
NPROCS=${NPROCS:-4}
