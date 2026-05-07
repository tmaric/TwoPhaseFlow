#!/bin/sh

# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
SOLVER=${SOLVER:-acousticHelmholtzFoam}

# ---------------------------------------------------------------------------
# Drive parameters
# ---------------------------------------------------------------------------
DRIVE_F=${DRIVE_F:-25250}
PISTON_U=${PISTON_U:-1.0}

# ---------------------------------------------------------------------------
# Fluid / acoustic medium properties
# ---------------------------------------------------------------------------
SOUND_SPEED=${SOUND_SPEED:-343}
# Default to a solid-like particle proxy for the ball study.
RHOL=${RHOL:-2500}
RHOG=${RHOG:-1.2}
CG=${CG:-343}
CL=${CL:-5000}

# ---------------------------------------------------------------------------
# PML parameters
# ---------------------------------------------------------------------------
PML_L=${PML_L:-0.008}
SIGMA_MAX=${SIGMA_MAX:-2000000}
PO=${PO:-3}
PML_MAX_Z=${PML_MAX_Z:-100}
PML_MIN_X=${PML_MIN_X:--100}
PML_MIN_Z=${PML_MIN_Z:--100}

# ---------------------------------------------------------------------------
# Acoustic / geometry parameters used to render the axisymmetric cross-section.
# ---------------------------------------------------------------------------
# Andrade 2019 wedge levitator geometry [m].
D=${D:-0.00764}
RT=${RT:-0.01}
# Use the original Andrade geometry with a larger reflector.
RR=${RR:-0.019}
BU=${BU:-0.01}
HS=${HS:-0.01}

# Wedge opening angle [deg].
WEDGE_DEG=${WEDGE_DEG:-1.0}

# Small nonzero radius used to avoid a degenerate axis surface during extrusion.
AXIS_RADIUS=${AXIS_RADIUS:-1.0e-5}
# Dummy thickness for the 2D cfMesh surface ribbon.
SURFACE_THICKNESS=${SURFACE_THICKNESS:-0.1}

# Drop definition [m].
DROP_RADIUS=${DROP_RADIUS:-0.001}
DROP_CENTER_Y=${DROP_CENTER_Y:-0.0040003}

# cfMesh controls [m].
CFMESH_MAX_CELL_SIZE=${CFMESH_MAX_CELL_SIZE:-3e-4}
CFMESH_DROP_CELL_SIZE=${CFMESH_DROP_CELL_SIZE:-2.0e-5}
# Interface-only semicircle refinement controls.
CFMESH_INTERFACE_SEGMENTS=${CFMESH_INTERFACE_SEGMENTS:-48}
