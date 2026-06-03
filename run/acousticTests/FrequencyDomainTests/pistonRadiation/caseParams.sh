#!/bin/sh

# Single source of truth for pistonRadiation setup.
# Values in this file are injected into:
# - constant/transportProperties
# - constant/pistonRadiation.geo
# - 0.orig/Pim
#
# Workflow:
# 1) Edit values here.
# 2) Run ./Allrun (or ./prepareCase to only regenerate inputs).

# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
# Frequency-domain acoustic solver.
SOLVER=acousticHelmholtzFoam

# ---------------------------------------------------------------------------
# Drive parameters
# ---------------------------------------------------------------------------
# Driven acoustic frequency [Hz].
DRIVE_F=10000
# Piston normal velocity amplitude used in 0.orig/Pim gradient BC [m/s].
PISTON_U=1

# ---------------------------------------------------------------------------
# Fluid / acoustic medium properties
# ---------------------------------------------------------------------------
# Liquid density [kg/m^3].
RHOL=998.3
# Gas density [kg/m^3].
RHOG=1.2
# Gas sound speed [m/s].
CG=343
# Liquid sound speed [m/s].
CL=1480

# ---------------------------------------------------------------------------
# PML parameters (spherical PML)
# ---------------------------------------------------------------------------
# Maximum damping strength sigmaMax [1/s].
# Tune per mesh/frequency setup.
SIGMA_MAX=100000 
# Polynomial profile order (po), typically 2-4.
PO=3
# Shifted-Laplacian strength beta (dimensionless).
# Used by freITBCFoam preconditioner only.
SHIFTED_LAPLACIAN_BETA=0.2
# Inner radius of PML region [m].
PML_RMIN=0.2
# PML thickness [m]. Outer radius is Rmin + L.
PML_L=0.08

# ---------------------------------------------------------------------------
# Mesh parameters (gmsh template)
# ---------------------------------------------------------------------------
# Number of cells per wavelength at DRIVE_F.
MESH_CELLS_PER_LAMBDA=${MESH_CELLS_PER_LAMBDA:-60}
# Transfinite progression ratio.
# 1.0 gives uniform spacing; >1.0 biases cell size growth along transfinite edges.
MESH_GR=1.0
# Wedge opening angle [deg] for axisymmetric approximation.
WEDGE_ANGLE_DEG=2
