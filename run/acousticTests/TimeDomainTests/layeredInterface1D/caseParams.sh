#!/bin/sh

# Solver (MPI default)
: "${SOLVER:=acousticWaveFoam}"

# Drive
: "${DRIVE_F:=20000}"
: "${PISTON_U:=0.01}"

# Medium 1 (x < X_INTERFACE)
: "${RHOG:=1.2}"
: "${CG:=343}"

# Medium 2 (x >= X_INTERFACE)
: "${RHOL:=1000}"
: "${CL:=1500}"

# Geometry [m]
: "${X_MIN:=0.0}"
: "${X_MAX:=0.35}"
: "${Y_HALF:=0.001}"
: "${Z_HALF:=0.001}"
: "${X_INTERFACE:=0.09}"

# Right PML [m]
: "${PML_L:=0.15}"
: "${PML_XMAX:=0.35}"
: "${SIGMA_MAX:=200000}"
: "${PO:=3}"


# Mesh
: "${NX:=2000}"
: "${NY:=1}"
: "${NZ:=1}"
