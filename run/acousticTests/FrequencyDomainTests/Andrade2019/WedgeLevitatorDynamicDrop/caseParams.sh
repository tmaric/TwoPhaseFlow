#!/bin/sh

# Single source of truth for WedgeLevitatorDropAlpha setup.
# Values in this file are injected into:
# - constant/transportProperties
# - constant/levitatorWedgeHex.geo
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
PISTON_U=2

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
SIGMA_MAX=5000
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

# Mesh controls.
# N is the number of cells per wavelength used to set gs = lambda/N.
MESH_N=100
# Gmsh point characteristic length (kept as in original setup).
MESH_LC=1.0
# Wedge half/total angle settings [deg].
ROTATE_HALF_DEG=0.5
WEDGE_DEG=1.0

# ---------------------------------------------------------------------------
# Initial drop setup (setAlphaField)
# ---------------------------------------------------------------------------
# Initial spherical drop radius [m].
DROP_RADIUS=0.0015
# Initial drop center vertical position [m] (y-coordinate).
DROP_CENTER_Y=0.004
