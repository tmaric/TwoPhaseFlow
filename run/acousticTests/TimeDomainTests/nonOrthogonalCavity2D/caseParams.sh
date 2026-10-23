#!/bin/sh

# Closed rigid square cavity with the exact standing mode
#   p = PA cos(MODE_M pi x/L) cos(MODE_N pi y/L) cos(omega t),
# on a mesh whose interior vertices are displaced smoothly,
#   x' = x + DISTORTION L/(2 pi) sin(pi x/L) sin(2 pi y/L),
#   y' = y + DISTORTION L/(2 pi) sin(pi y/L) sin(2 pi x/L),
# which keeps the walls straight and makes the faces non-orthogonal.
# No PML, no interface, no damping. Every entry can be overridden.

: "${SOLVER:=acousticWaveFoam}"
: "${RHOG:=1.2}"
: "${CG:=343}"
: "${CAVITY_L:=0.343}"
: "${MODE_M:=2}"
: "${MODE_N:=1}"
: "${PA:=1}"

# Mesh: CELLS per side; DISTORTION 0 gives the orthogonal mesh.
: "${CELLS:=40}"
: "${DISTORTION:=0}"
# WALL_LAYERS cell layers next to each wall stay orthogonal; the distortion
# then ramps up smoothly over TRANSITION_FRACTION L (as in the frequency
# case homogeneousPlaneWave2D). 0 and 0 distort up to the walls.
: "${WALL_LAYERS:=0}"
: "${TRANSITION_FRACTION:=0}"
# Laplacian and surface-normal gradient: orthogonal | corrected. With
# corrected, the non-orthogonal correction is an explicit source: lagged by
# one step in the monolithic update, iterated NONORTH_CORRECTORS times within
# the step in the segregated update.
: "${LAPLACIAN:=orthogonal}"
: "${NONORTH_CORRECTORS:=0}"
# Monolithic update with LAPLACIAN=uncorrected: none | explicit (the
# correction of the corrected scheme, applied by the solver) | symmetric
# (its symmetrized form, which keeps the operator self-adjoint).
: "${NONORTH_MODE:=none}"
# Monolithic: stop the correctors at this relative pressure change (0: run
# all NONORTH_CORRECTORS).
: "${NONORTH_TOLERANCE:=0}"

# Time: STEPS_PER_PERIOD of the continuous mode, N_PERIODS periods,
# PROBE_SAMPLES_PER_PERIOD samples per period.
: "${STEPS_PER_PERIOD:=160}"
: "${N_PERIODS:=20}"
: "${PROBE_SAMPLES_PER_PERIOD:=8}"

# blockCoupled (monolithic) | thetaNewmark (segregated) | fullyImplicit
: "${TIME_METHOD:=blockCoupled}"
: "${THETA:=0.5}"
: "${BLOCK_SOLVER:=direct}"
