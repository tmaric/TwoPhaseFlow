#!/bin/sh

# Single source of truth for WedgeLevitatorDropHole setup.
# This case follows WedgeLevitatorDropAlpha except that the drop is meshed as
# a sound-hard boundary patch named dropWall instead of an alpha.water region.

# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
# Frequency-domain acoustic solver:
# - acousticHelmholtzFoam: MPI-capable
# - acousticHelmholtzSerialFoam: single-rank reference solver
SOLVER=acousticHelmholtzFoam

# ---------------------------------------------------------------------------
# Drive parameters
# ---------------------------------------------------------------------------
# Driven acoustic frequency [Hz].
DRIVE_F=25250
# Piston normal velocity amplitude for transducer gradient BC [m/s].
PISTON_U=1.0

# ---------------------------------------------------------------------------
# Fluid / acoustic medium properties
# ---------------------------------------------------------------------------
# Gas sound speed used by geometry and wave-number definitions [m/s].
SOUND_SPEED=343
# Liquid density [kg/m^3]. Kept for solver consistency with DropAlpha.
RHOL=998.3
# Gas density [kg/m^3].
RHOG=1.2
# Gas sound speed [m/s].
CG=343
# Liquid sound speed [m/s]. Kept for solver consistency with DropAlpha.
CL=1480

# ---------------------------------------------------------------------------
# PML parameters (rectangle PML)
# ---------------------------------------------------------------------------
# PML thickness [m].
PML_L=0.008
# Max damping strength sigmaMax [1/s]. tLF:1E7 aH(S)F:2E6
SIGMA_MAX=2000000
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
# Levitator geometry dimensions [m], matching DropAlpha.
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

# Wedge angle [deg], matching DropAlpha.
WEDGE_DEG=1.0

# Small nonzero radius used to avoid a degenerate axis surface before extrusion.
AXIS_RADIUS=${AXIS_RADIUS:-1.0e-5}
# Dummy thickness for the 2D cfMesh surface ribbon.
SURFACE_THICKNESS=${SURFACE_THICKNESS:-0.1}

# Sound-hard drop boundary [m].
DROP_RADIUS=0.0001
DROP_HORIZONTAL_LONG_AXIS=0.0001
DROP_CENTER_X=1e-8
DROP_CENTER_Y=0.003550003
DROP_CENTER_Z=1e-8

# cfMesh controls [m]. These replace DropAlpha's gmsh/setAlphaField interface
# handling and control the hard dropWall boundary resolution.
CFMESH_MAX_CELL_SIZE=${CFMESH_MAX_CELL_SIZE:-1e-4}
CFMESH_DROP_CELL_SIZE=${CFMESH_DROP_CELL_SIZE:-5e-5}
CFMESH_DROP_SEGMENTS=${CFMESH_DROP_SEGMENTS:-128}
