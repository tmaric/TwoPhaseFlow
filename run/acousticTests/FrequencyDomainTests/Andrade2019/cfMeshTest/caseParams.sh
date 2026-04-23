#!/bin/sh

# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
SOLVER=acousticHelmholtzFoam

# ---------------------------------------------------------------------------
# Drive parameters
# ---------------------------------------------------------------------------
DRIVE_F=25250
PISTON_U=1.0

# ---------------------------------------------------------------------------
# Fluid / acoustic medium properties
# ---------------------------------------------------------------------------
SOUND_SPEED=343
RHOL=998.3
RHOG=1.2
CG=343
CL=1480

# ---------------------------------------------------------------------------
# PML parameters
# ---------------------------------------------------------------------------
PML_L=0.008
SIGMA_MAX=2000000
PO=3
PML_MAX_Z=100
PML_MIN_X=-100
PML_MIN_Z=-100

# ---------------------------------------------------------------------------
# Acoustic / geometry parameters used to render the axisymmetric cross-section.
# ---------------------------------------------------------------------------
# Andrade 2019 wedge levitator geometry [m].
D=0.00764
RT=0.01
# Use the original Andrade geometry with a larger reflector.
RR=0.019
BU=0.01
HS=0.01

# Wedge opening angle [deg].
WEDGE_DEG=1.0

# Small nonzero radius used to avoid a degenerate axis surface during extrusion.
AXIS_RADIUS=1.0e-5
# Dummy thickness for the 2D cfMesh surface ribbon.
SURFACE_THICKNESS=0.1

# Drop definition [m].
DROP_RADIUS=0.001
DROP_CENTER_Y=0.0040003

# cfMesh controls [m].
CFMESH_MAX_CELL_SIZE=1.0e-4
CFMESH_DROP_CELL_SIZE=3.0e-5
# Interface-only semicircle refinement controls.
CFMESH_INTERFACE_SEGMENTS=24
