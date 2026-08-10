#!/bin/sh

# Solver (MPI default)
SOLVER=acousticHelmholtzFoam
CASE_MODE=${CASE_MODE:-pml}

# Drive
DRIVE_F=${DRIVE_F:-20000}
PISTON_U=${PISTON_U:-0.01}

# Medium 1 (x < X_INTERFACE)
RHOG=${RHOG:-1000}
CG=${CG:-1500}

# Medium 2 (x >= X_INTERFACE)
RHOL=${RHOL:-1.2}
CL=${CL:-343}

# Geometry [m]
X_MIN=${X_MIN:-0.0}
X_MAX=${X_MAX:-0.35}
Y_HALF=0.001
Z_HALF=0.001
X_INTERFACE=${X_INTERFACE:-0.09}

# Right PML [m]
PML_L=${PML_L:-0.15}
PML_XMAX=${PML_XMAX:-$X_MAX}
SIGMA_MAX=${SIGMA_MAX:-500000}
PO=${PO:-3}

# Mesh
NX=${NX:-8000}
NY=1
NZ=1
