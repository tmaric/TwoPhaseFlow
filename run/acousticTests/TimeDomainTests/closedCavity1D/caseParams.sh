#!/bin/sh

# Closed rigid one-dimensional cavity with the exact standing mode
#   p = PA cos(k x) cos(omega t),  k = MODE_N pi / CAVITY_L,  omega = CG k.
# No PML, no material interface and no added damping. Every entry can be
# overridden from the environment.

: "${SOLVER:=acousticWaveFoam}"

# Single gas medium; alpha.water is absent, so alpha = 0 everywhere.
: "${RHOG:=1.2}"
: "${CG:=343}"

# Cavity and mode
: "${CAVITY_L:=0.343}"
: "${MODE_N:=2}"
: "${PA:=1}"

# Resolution. CELLS_PER_WAVELENGTH*MODE_N/2 must be an integer.
: "${CELLS_PER_WAVELENGTH:=40}"
: "${STEPS_PER_PERIOD:=20}"
: "${N_PERIODS:=100}"
# Probe samples per period; STEPS_PER_PERIOD must be a multiple of it.
# The backward (fullyImplicit) discrete energy needs every step.
: "${PROBE_SAMPLES_PER_PERIOD:=20}"

# Time integration: blockCoupled | thetaNewmark | fullyImplicit
: "${TIME_METHOD:=blockCoupled}"
: "${THETA:=0.5}"
# blockCoupled only: direct | iterative
: "${BLOCK_SOLVER:=direct}"

# Radiation averages: window length in periods of the continuous mode (the
# solver rounds it to an integer number of steps). RADIATION_CHECK=1 runs
# postprocess_radiation.py, which needs every step probed.
: "${RADIATION_PERIODS:=1}"
: "${RADIATION_CHECK:=0}"
