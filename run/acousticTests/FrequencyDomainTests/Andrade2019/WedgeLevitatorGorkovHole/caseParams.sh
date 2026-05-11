#!/bin/sh

# Single source of truth for WedgeLevitatorGorkovHole setup.
# This case follows the GorkovAlpha postprocessing workflow, but uses the
# DropHole mesh treatment: the drop is a sound-hard boundary patch named
# dropWall instead of an alpha.water region.

# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
# Frequency-domain acoustic solver:
# - acousticHelmholtzFoam: MPI-capable
# - acousticHelmholtzSerialFoam: single-rank reference solver
SOLVER=${SOLVER:-acousticHelmholtzFoam}

# ---------------------------------------------------------------------------
# Drive parameters
# ---------------------------------------------------------------------------
# Driven acoustic frequency [Hz].
DRIVE_F=${DRIVE_F:-25250}
# Piston normal velocity amplitude for transducer gradient BC [m/s].
PISTON_U=${PISTON_U:-1.0}

# ---------------------------------------------------------------------------
# Fluid / acoustic medium properties
# ---------------------------------------------------------------------------
# Gas sound speed used by geometry and wave-number definitions [m/s].
SOUND_SPEED=${SOUND_SPEED:-343}
# Particle density [kg/m^3].
RHOL=${RHOL:-2500}
# Gas density [kg/m^3].
RHOG=${RHOG:-1.2}
# Gas sound speed [m/s].
CG=${CG:-343}
# Particle sound speed [m/s].
CL=${CL:-5000}

# ---------------------------------------------------------------------------
# PML parameters (rectangle PML)
# ---------------------------------------------------------------------------
# PML thickness [m].
PML_L=${PML_L:-0.008}
# Max damping strength sigmaMax [1/s]. tLF:1E7 aH(S)F:2E6
SIGMA_MAX=${SIGMA_MAX:-2000000}
# Polynomial profile order (po), typically 2-4.
PO=${PO:-3}
# Rectangle bounds for PML logic in transportProperties.
# Derived automatically in prepareCase:
#   PML_MAX_X = BU + RR
#   PML_MAX_Y = HS + D
#   PML_MIN_Y = -HS
# Keep only fixed limits here.
PML_MAX_Z=${PML_MAX_Z:-100}
PML_MIN_X=${PML_MIN_X:--100}
PML_MIN_Z=${PML_MIN_Z:--100}

# ---------------------------------------------------------------------------
# Geometry / mesh parameters
# ---------------------------------------------------------------------------
# Levitator geometry dimensions [m], matching DropAlpha.
# D: gap between transducer plane and reflector plane.
D=${D:-0.00764}
# RT: transducer radius.
RT=${RT:-0.01}
# RR: reflector radius.
RR=${RR:-0.019}
# BU: radial air buffer beyond reflector radius (outer side extension).
BU=${BU:-0.01}
# HS: vertical stand height above/below the active gap.
HS=${HS:-0.01}

# Wedge angle [deg], matching DropAlpha.
WEDGE_DEG=${WEDGE_DEG:-1.0}

# Small nonzero radius used to avoid a degenerate axis surface before extrusion.
AXIS_RADIUS=${AXIS_RADIUS:-1.0e-5}
# Dummy thickness for the 2D cfMesh surface ribbon.
SURFACE_THICKNESS=${SURFACE_THICKNESS:-0.1}

# Sound-hard drop boundary [m].
DROP_RADIUS=${DROP_RADIUS:-0.00002}
DROP_CENTER_Y=${DROP_CENTER_Y:-0.0040003}

# cfMesh controls [m]. These replace DropAlpha's gmsh/setAlphaField interface
# handling and control the hard dropWall boundary resolution.
CFMESH_MAX_CELL_SIZE=${CFMESH_MAX_CELL_SIZE:-1e-4}
CFMESH_DROP_CELL_SIZE=${CFMESH_DROP_CELL_SIZE:-5e-6}
CFMESH_DROP_SEGMENTS=${CFMESH_DROP_SEGMENTS:-256}

# Optional per-run overrides used by parameter sweeps. This keeps copied sweep
# cases self-contained when they are rerun or inspected later.
if [ -f ./caseOverrides.sh ]; then
    . ./caseOverrides.sh
fi
