#!/bin/sh

# Solver (MPI default)
SOLVER=acousticHelmholtzFoam

# Drive
DRIVE_F=10000
PISTON_U=0.01

# Medium 1 (x < X_INTERFACE)
RHOG=1.2
CG=343

# Medium 2 (x >= X_INTERFACE)
RHOL=20
CL=343

# Geometry [m]
X_MIN=0.0
X_MAX=0.30
Y_HALF=0.001
Z_HALF=0.001
X_INTERFACE=0.12

# Right PML [m]
PML_L=0.06
PML_XMAX=0.30
SIGMA_MAX=1000
PO=3

# Mesh
NX=600
NY=1
NZ=1
