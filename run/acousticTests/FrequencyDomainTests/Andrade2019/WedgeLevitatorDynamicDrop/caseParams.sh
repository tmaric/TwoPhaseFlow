#!/bin/sh

# Single source of truth for WedgeLevitatorDropAlpha setup.
# Values in this file are injected into:
# - constant/transportProperties
# - 0.orig/Pim
#
# Workflow:
# 1) Edit values here.
# 2) Run ./Allrun (or ./prepareCase to only regenerate inputs).

# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
# MPI-capable coupled solver.
SOLVER=interFALFlow

# ---------------------------------------------------------------------------
# Drive parameters
# ---------------------------------------------------------------------------
# Driven acoustic frequency [Hz].
DRIVE_F=25250
# Piston normal velocity amplitude for transducer gradient BC [m/s].
PISTON_U=1

# ---------------------------------------------------------------------------
# Fluid / acoustic medium properties
# ---------------------------------------------------------------------------
# Gas sound speed used by geometry and wave-number definitions [m/s].
SOUND_SPEED=343
# Liquid density [kg/m^3].
RHOL=998.3
# Gas density [kg/m^3].
RHOG=1.2
# Liquid kinematic viscosity [m^2/s].
NUL=1.004e-6
# Gas kinematic viscosity [m^2/s].
NUG=1.51e-5
# Gas sound speed [m/s].
CG=343
# Liquid sound speed [m/s].
CL=1480
# ---------------------------------------------------------------------------
# PML parameters (rectangle PML)
# ---------------------------------------------------------------------------
# PML thickness [m].
PML_L=0.008
# Directional maximum damping strengths [1/s].
# The wedge uses rectangular PML regions in x and y, but not in z.
SIGMA_MAX=${SIGMA_MAX:-2000000}
SIGMA_MAX_X=${SIGMA_MAX_X:-${SIGMA_MAX}}
SIGMA_MAX_Y=${SIGMA_MAX_Y:-${SIGMA_MAX}}
SIGMA_MAX_Z=${SIGMA_MAX_Z:-0}
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
# Geometry / mesh parameters
# ---------------------------------------------------------------------------
# Levitator geometry dimensions [m].
# D: gap between transducer plane and reflector plane.
D=0.00764
# RT: transducer radius.
RT=0.01
# RR: reflector radius.
RR=0.019
# BU: radial air buffer beyond reflector radius (outer side extension).
BU=0.01
# HS: vertical stand height above/below the active gap.
HS=0.01

# Wedge total angle [deg].
WEDGE_DEG=1.0

# cfMesh controls [m].
AXIS_RADIUS=${AXIS_RADIUS:-1.0e-5}
SURFACE_THICKNESS=${SURFACE_THICKNESS:-0.1}
CFMESH_MAX_CELL_SIZE=${CFMESH_MAX_CELL_SIZE:-3e-5}
CFMESH_INTERFACE_REFINEMENT=${CFMESH_INTERFACE_REFINEMENT:-0}
CFMESH_INTERFACE_CELL_SIZE=${CFMESH_INTERFACE_CELL_SIZE:-3e-5}
CFMESH_INTERFACE_BOX_SCALE=${CFMESH_INTERFACE_BOX_SCALE:-2}

# ---------------------------------------------------------------------------
# Initial drop setup (setAlphaField)
# ---------------------------------------------------------------------------
# Volume-equivalent spherical drop radius [m].
DROP_RADIUS=0.001
# Horizontal long semi-axis [m]. Set equal to DROP_RADIUS for a sphere.
DROP_HORIZONTAL_LONG_AXIS=0.001
# Small x/z offsets keep the center away from the wedge singular line.
DROP_CENTER_X=1e-8
# Initial drop center vertical position [m] (y-coordinate).
DROP_CENTER_Y=0.004
DROP_CENTER_Z=1e-8
