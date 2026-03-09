#!/bin/sh

# Single source of truth for WedgeLevitatorHeightSweep setup.
# Values in this file are injected into:
# - constant/transportProperties
# - constant/levitatorWedgeHex.geo
# - 0.orig/Pim
#
# Workflow:
# 1) Edit values here.
# 2) Run ./Allrun (or ./prepareCase to only regenerate inputs).

# ---------------------------------------------------------------------------
# Run mode
# ---------------------------------------------------------------------------
# Solver mode used by Allrun:
# - serial: single-rank acoustic solve
# - mpi:    decomposed/multi-rank acoustic solve
: "${RUN_MODE:=serial}"
SERIAL_SOLVER=acousticHelmholtzFoam
MPI_SOLVER=freBCMFoam

# ---------------------------------------------------------------------------
# Drive parameters
# ---------------------------------------------------------------------------
# Driven acoustic frequency [Hz].
DRIVE_F=25250
# Piston normal velocity amplitude for transducer gradient BC [m/s].
PISTON_U=1.0

# Reflector height normalized by wavelength: D/lambda.
# The geometric gap D is computed in prepareCase as:
#   D = HEIGHT_FAC * (SOUND_SPEED/DRIVE_F)
: "${HEIGHT_FAC:=1.0}"

# ---------------------------------------------------------------------------
# Fluid / acoustic medium properties
# ---------------------------------------------------------------------------
# Gas sound speed used by geometry and wave-number definitions [m/s].
SOUND_SPEED=343
# Liquid density [kg/m^3].
RHOL=998.3
# Gas density [kg/m^3].
RHOG=1.2
# Gas sound speed [m/s].
CG=343
# Liquid sound speed [m/s].
CL=1480

# ---------------------------------------------------------------------------
# PML parameters (rectangle PML)
# ---------------------------------------------------------------------------
# PML thickness [m].
PML_L=0.008
# Max damping strength sigmaMax [1/s].
SIGMA_MAX=3000
# Polynomial profile order (po), typically 2-4.
PO=3
# Rectangle bounds for PML logic in transportProperties.
# Derived automatically in prepareCase:
#   PML_MAX_X = BU + RR
#   PML_MAX_Y = HS + D
#   PML_MIN_Y = -HS
# Keep only fixed limits here.
PML_MAX_Z=100
PML_MIN_X=-100
PML_MIN_Z=-100

# ---------------------------------------------------------------------------
# Geometry / mesh parameters (gmsh template)
# ---------------------------------------------------------------------------
# Levitator geometry dimensions [m].
# RT: transducer radius.
RT=0.01
# RR: reflector radius.
RR=0.019
# BU: radial air buffer beyond reflector radius (outer side extension).
BU=0.01
# HS: vertical stand height above/below the active gap.
HS=0.01

# Mesh controls.
# N is the number of cells per wavelength used to set gs = lambda/N.
MESH_N=200
# Gmsh point characteristic length (kept as in original setup).
MESH_LC=1.0
# Wedge half/total angle settings [deg].
ROTATE_HALF_DEG=0.5
WEDGE_DEG=1.0

# ---------------------------------------------------------------------------
# Height sweep controls
# ---------------------------------------------------------------------------
# Sweep range in normalized height D/lambda.
: "${SWEEP_FAC_MIN:=0.4}"
: "${SWEEP_FAC_MAX:=1.6}"
: "${SWEEP_POINTS:=100}"
