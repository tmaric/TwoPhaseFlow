# pistonPulse1D

A piston at `x = 0` emits a smooth velocity pulse that crosses a planar
gas|liquid interface and leaves the physical domain through a PML, solved
with `acousticWaveFoam`. The case verifies time-dependent boundary data, the
interface and the PML together against an exact solution, with zero added
damping.

## Configuration

- Geometry of the layered case: medium 1 on `0 < x < X_INTERFACE`
  (0.09 m), medium 2 up to `LENGTH` (0.35 m), and a cubic Cartesian PML of
  thickness `PML_THICKNESS` (0.15 m) and `SIGMA_MAX` (5e5 1/s) at the end of
  medium 2, terminated by a sound-hard wall.
- `CONFIG = gasFirst` puts the piston in air and the PML in water;
  `liquidFirst` reverses the media.
- Piston velocity `U0 sin(w0 t) sin^4(pi t/Tp)` for `0 <= t <= Tp`,
  `Tp = N_CYCLES/F0` (3 cycles at 20 kHz), zero afterwards. It enters as
  the boundary gradient `dp/dn = rho1 dU_p/dt` along the outward normal, as
  an OpenFOAM `expression` Function1 in time.
- `NX` must place the interface on a face (a multiple of 35 for the
  defaults); `STEPS_PER_PERIOD` counts time steps per period `1/F0`.

## Exact solution

`pulseReference.py` evaluates the finite sum of travelling pulses for an
ideal PML: with `R = (Z2 - Z1)/(Z2 + Z1)` the pulse bounces between the
interface and the piston, which reflects like a rigid wall, and each passage
transmits a pulse into medium 2. For air and water `R = 0.99945`.

## Run

Load OpenFOAM-v2606 and `scripts/bashrc`, build `acousticWaveFoam`, then

```bash
./Allrun
```

Every entry of `caseParams.sh` can be overridden, e.g.
`CONFIG=liquidFirst NX=1400 STEPS_PER_PERIOD=700 ./Allrun`.

## Post-processing

`postprocess_pulse.py` compares the probed pressure and velocity with the
exact solution over the physical domain, with energy weights (`V beta` for
pressure, `V rho` for velocity), for the whole domain and per medium. It
also reports the largest pressure in the physical part of medium 2 while the
exact pressure vanishes there, relative to the transmitted amplitude: the
reflection from the PML plus numerical tails. Results:
`postProcessing/pulseVerification/`.

## Studies

```bash
./runStudy A     # temporal refinement, monolithic and backward (10 runs)
./runStudy B     # spatial refinement, piston in gas (4 runs)
./runStudy C     # PML strength and thickness (8 runs)
./runStudy D     # spatial refinement, piston in liquid (4 runs)
./runStudy       # all
```

`summarize_pulse_study.py` writes `studies/study*.csv`.
