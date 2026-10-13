#!/bin/sh

# Gaussian pressure pulse in a square (DIM=2) or cube (DIM=3) of gas,
# surrounded by Cartesian PMLs on all sides, sound-hard outer boundary.
# Every entry can be overridden from the environment.

: "${SOLVER:=acousticWaveFoam}"

# Media; the drop (DROP_RADIUS > 0) is liquid.
: "${RHOG:=1.2}"
: "${CG:=343}"
: "${RHOL:=1000}"
: "${CL:=1500}"

# Geometry [m]: physical region |x_i| < HALF_WIDTH, PML of PML_THICKNESS
# outside it. ENLARGED=1 removes the PML and extends the domain to
# ENLARGED_HALF_WIDTH (reference run); (ENLARGED_HALF_WIDTH - HALF_WIDTH -
# PML_THICKNESS)/CELL_SIZE must be an integer so that the cell centres of
# both grids coincide.
: "${DIM:=2}"
: "${HALF_WIDTH:=0.05}"
: "${PML_THICKNESS:=0.025}"
: "${SIGMA_MAX:=5.05e5}"
: "${PML_ORDER:=3}"
: "${ENLARGED:=0}"
: "${ENLARGED_HALF_WIDTH:=0.115}"
: "${CELL_SIZE:=1e-3}"
: "${DROP_RADIUS:=0}"

# Initial Gaussian pulse P0 exp(-|x - x0|^2/(2 s^2)), zero rate
: "${P0:=1}"
: "${PULSE_WIDTH:=4e-3}"
: "${SOURCE:=0.015 0.01 0}"

# Time: COURANT = cg dt / h; samples every SAMPLE_INTERVAL [s]
: "${COURANT:=0.5}"
: "${END_TIME:=4e-4}"
: "${SAMPLE_INTERVAL:=1e-5}"
# 1: also probe p, psi (v) and w (u) at every cell, PML included, for the
# long-time maxima; 0 for fine accuracy runs.
: "${DOMAIN_PROBES:=1}"

# Time integration: blockCoupled (monolithic) | thetaNewmark (segregated)
: "${TIME_METHOD:=blockCoupled}"
: "${THETA:=0.5}"
: "${BLOCK_SOLVER:=direct}"
