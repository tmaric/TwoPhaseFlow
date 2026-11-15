#!/bin/sh

# Canonical rigid-sphere verification in a one-dimensional standing wave.
# The analytical amplitude is imposed only at the transducer boundary; the
# pressure around the sound-hard sphere is obtained from the numerical solve.

DRIVE_F=${DRIVE_F:-25250}
PRESSURE_AMPLITUDE=${PRESSURE_AMPLITUDE:-1.0}

RHOG=${RHOG:-1.2}
CG=${CG:-343}
RHOL=${RHOL:-2500}
CL=${CL:-5000}

# Equal transducer and reflector radii give a cylindrical one-dimensional
# reference cavity. The side wall and reflector are sound hard.
H=${H:-0.00764}
CAVITY_RADIUS=${CAVITY_RADIUS:-0.01}
WEDGE_DEG=${WEDGE_DEG:-1.0}
AXIS_RADIUS=${AXIS_RADIUS:-1.0e-8}
SURFACE_THICKNESS=${SURFACE_THICKNESS:-0.1}

DROP_RADIUS=${DROP_RADIUS:-100e-6}
DROP_CENTER_Y=${DROP_CENTER_Y:-0.0040003}

CFMESH_MAX_CELL_SIZE=${CFMESH_MAX_CELL_SIZE:-6e-5}
CFMESH_DROP_CELL_SIZE=${CFMESH_DROP_CELL_SIZE:-2e-5}
CFMESH_DROP_SEGMENTS=${CFMESH_DROP_SEGMENTS:-64}

# The solver requires a sigma field. Zero damping makes this a closed-cavity
# verification without an active PML.
PML_L=${PML_L:-0.001}
SIGMA_MAX=${SIGMA_MAX:-0}
PO=${PO:-3}

if [ -f ./caseOverrides.sh ]; then
    . ./caseOverrides.sh
fi
