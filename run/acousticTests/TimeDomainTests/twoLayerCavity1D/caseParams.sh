#!/bin/sh

# Closed rigid one-dimensional cavity with a gas layer on [0, x_I] and a
# liquid layer on [x_I, L]; no PML, no forcing and no added damping. The
# initial pressure is an exact eigenmode of the two-layer problem. Every
# entry can be overridden from the environment.

: "${SOLVER:=acousticWaveFoam}"

# Media
: "${RHOG:=1.2}"
: "${CG:=343}"
: "${RHOL:=1000}"
: "${CL:=1500}"

# Geometry in units u = L_GAS/GAS_UNITS. The defaults place the second gas
# half-wave resonance and the first liquid quarter-wave resonance near 1 kHz,
# so the two modes there hybridize and carry energy in both layers.
: "${L_GAS:=0.343}"
: "${GAS_UNITS:=32}"
: "${LIQUID_UNITS:=35}"
# Interface shift into the liquid as a fraction of the cell size h;
# 0 places it on a face, 0 < s < 1 inside a cell.
: "${INTERFACE_OFFSET:=0}"

# Mode: the exact eigenmode nearest to F_TARGET from below or above.
: "${F_TARGET:=1000}"
: "${MODE:=lower}"
: "${PA:=1}"

# Resolution: cells per unit u (the gas wavelength is 32 u by default).
: "${CELLS_PER_UNIT:=2}"
: "${STEPS_PER_PERIOD:=20}"
: "${N_PERIODS:=100}"
# Probe samples per period; STEPS_PER_PERIOD must be a multiple of it.
: "${PROBE_SAMPLES_PER_PERIOD:=20}"

# Time integration: blockCoupled | thetaNewmark | fullyImplicit
: "${TIME_METHOD:=blockCoupled}"
: "${THETA:=0.5}"
# blockCoupled only: direct | iterative
: "${BLOCK_SOLVER:=direct}"

# 1: serial. 2: two ranks split at the interface; with an offset the
# interface cell is the first cell of rank 1, next to the processor face.
: "${NPROCS:=1}"

# Radiation averages: window length in periods of the mode (the solver
# rounds it to an integer number of steps). RADIATION_CHECK=1 runs
# postprocess_radiation.py (serial runs with the interface on a face).
: "${RADIATION_PERIODS:=1}"
: "${RADIATION_CHECK:=0}"
