#!/bin/sh

# Single source of truth for the empty-levitator setup used by DropHole.
# This case uses the same gmsh/gmshToFoam outer levitator mesh workflow as
# DropAlpha, but omits the drop/hole. A function object writes the Gorkov
# potential field for a small probe particle in the empty acoustic field.

# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
SOLVER=${SOLVER:-acousticHelmholtzFoam}

# ---------------------------------------------------------------------------
# Drive parameters
# ---------------------------------------------------------------------------
DRIVE_F=${DRIVE_F:-100000}
PISTON_U=${PISTON_U:-1.0}

# ---------------------------------------------------------------------------
# Host fluid and probe particle properties
# ---------------------------------------------------------------------------
SOUND_SPEED=${SOUND_SPEED:-343}
RHOG=${RHOG:-1.2}
CG=${CG:-343}

# Probe particle properties used only in the Gorkov-potential function object.
RHOL=${RHOL:-2500}
CL=${CL:-5000}
PARTICLE_RADIUS=${PARTICLE_RADIUS:-100e-6}
PARTICLE_CENTER_Y=${PARTICLE_CENTER_Y:-0.0040003}

# ---------------------------------------------------------------------------
# PML parameters
# ---------------------------------------------------------------------------
PML_L=${PML_L:-0.008}
SIGMA_MAX=${SIGMA_MAX:-1000000}
SIGMA_MAX_X=${SIGMA_MAX_X:-${SIGMA_MAX}}
SIGMA_MAX_Y=${SIGMA_MAX_Y:-${SIGMA_MAX}}
SIGMA_MAX_Z=${SIGMA_MAX_Z:-0}
PO=${PO:-3}
PML_MAX_Z=${PML_MAX_Z:-100}
PML_MIN_X=${PML_MIN_X:--100}
PML_MIN_Z=${PML_MIN_Z:--100}

# ---------------------------------------------------------------------------
# Geometry / mesh parameters
# ---------------------------------------------------------------------------
D=${D:-0.00764}
RT=${RT:-0.01}
RR=${RR:-0.019}
BU=${BU:-0.01}
HS=${HS:-0.01}
WEDGE_DEG=${WEDGE_DEG:-1.0}

AXIS_RADIUS=${AXIS_RADIUS:-1.0e-5}
MESH_N=${MESH_N:-50}
MESH_LC=${MESH_LC:-1.0}
ROTATE_HALF_DEG=${ROTATE_HALF_DEG:-0.5}

if [ -f ./caseOverrides.sh ]; then
    . ./caseOverrides.sh
fi
