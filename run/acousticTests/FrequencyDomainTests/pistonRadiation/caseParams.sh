#!/bin/sh

# Single source of truth for pistonRadiation setup.
# Values in this file are injected into:
# - constant/transportProperties
# - system/blockMeshDict
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
# Exact piston radius [m]. The block boundary at this radius prevents the
# driven patch geometry from changing between mesh-convergence levels.
PISTON_RADIUS=0.1

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
# PML parameters (rectangular PML)
# ---------------------------------------------------------------------------
# Maximum damping strength sigmaMax [1/s].
# Tune per mesh/frequency setup.
SIGMA_MAX=100000 
# Polynomial profile order (po), typically 2-4.
PO=3
# Shifted-Laplacian strength beta (dimensionless).
# Used by freITBCFoam preconditioner only.
SHIFTED_LAPLACIAN_BETA=0.2
# Physical-domain extent in the radial and axial directions [m].
PML_RMIN=0.2
# Outer-domain extent in the radial and axial directions [m].
PML_RMAX=0.28

# ---------------------------------------------------------------------------
# Mesh parameters (blockMesh template)
# ---------------------------------------------------------------------------
# Number of cells per wavelength at DRIVE_F.
MESH_CELLS_PER_LAMBDA=${MESH_CELLS_PER_LAMBDA:-80}
# Overall block grading ratio in the axial and radial directions.
# 1.0 gives uniform spacing.
MESH_GR=1.0
# Wedge opening angle [deg] for axisymmetric approximation.
WEDGE_ANGLE_DEG=2
