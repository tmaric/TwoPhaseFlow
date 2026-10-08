#!/bin/sh

# Piston pulse through a gas|liquid interface into a PML, one dimension.
# Medium 1 (at the piston) fills 0 < x < X_INTERFACE, medium 2 the rest;
# the PML lies in medium 2 at the end of the domain. Every entry can be
# overridden from the environment.

: "${SOLVER:=acousticWaveFoam}"

# Media
: "${RHOG:=1.2}"
: "${CG:=343}"
: "${RHOL:=1000}"
: "${CL:=1500}"
# gasFirst: piston in gas, PML in liquid. liquidFirst: the reverse.
: "${CONFIG:=gasFirst}"

# Geometry [m]; NX must place X_INTERFACE on a face.
: "${LENGTH:=0.35}"
: "${X_INTERFACE:=0.09}"
: "${PML_THICKNESS:=0.15}"
: "${SIGMA_MAX:=5e5}"
: "${PML_ORDER:=3}"
: "${NX:=2800}"

# Piston pulse U0 sin(w0 t) sin^4(pi t/Tp), Tp = N_CYCLES/F0
: "${F0:=20000}"
: "${N_CYCLES:=3}"
: "${U0:=0.01}"

# Time: STEPS_PER_PERIOD per period 1/F0; observation interval [s]
: "${STEPS_PER_PERIOD:=1400}"
: "${T_OBS:=1.6e-3}"
# Probe samples per period; STEPS_PER_PERIOD must be a multiple of it.
: "${PROBE_SAMPLES_PER_PERIOD:=10}"

# Time integration: blockCoupled | thetaNewmark | fullyImplicit
: "${TIME_METHOD:=blockCoupled}"
: "${THETA:=0.5}"
: "${BLOCK_SOLVER:=direct}"
