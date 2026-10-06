# closedCavity1D

Free standing mode in a closed rigid one-dimensional cavity, solved with
`acousticWaveFoam`. The case has no PML, no material interface, no boundary
forcing and no added damping, so it isolates the time integration.

## Configuration

- Air (`RHOG`, `CG`); `alpha.water` is absent, so `alpha = 0` everywhere.
- Rigid walls (`zeroGradient` pressure) at `x = 0` and `x = CAVITY_L`.
- Initial state `p = PA cos(k x)`, `pDot = 0`, `U = 0` with
  `k = MODE_N pi / CAVITY_L`, which matches the exact mode
  `p = PA cos(k x) cos(omega t)` at `t = 0`.
- Uniform mesh with `CELLS_PER_WAVELENGTH` cells per wavelength and
  `STEPS_PER_PERIOD` time steps per period; `sigma` is zero everywhere.

The defaults (`CAVITY_L = 0.343 m`, `MODE_N = 2`) give a 1 kHz mode with the
wavelength equal to the cavity length.

## Run

Load OpenFOAM-v2606 and `scripts/bashrc`, build `acousticWaveFoam`, then

```bash
./Allrun
```

All entries of `caseParams.sh` can be overridden from the environment, e.g.

```bash
TIME_METHOD=fullyImplicit STEPS_PER_PERIOD=40 ./Allrun
```

`Allrun` runs serially: `blockMesh`, `setExprFields`, the solver, and
`postprocess_cavity.py`.

## Post-processing

Pressure, `pDot` and `U` are probed at every cell centre every
`STEPS_PER_PERIOD/PROBE_SAMPLES_PER_PERIOD` steps. On this mesh `cos(k x_P)`
is an exact eigenvector of the orthogonal pressure operator with frequency
`omega_h = (2 c / h) sin(k h / 2)`. The script projects the pressure onto this
mode and compares the amplitude with

- the exact amplitude `cos(omega t)`,
- the semi-discrete amplitude `cos(omega_h t)`,
- the fully discrete amplitude of the selected integrator
  (theta rule: `Re(g^n)`, `g = (1 + (1 - theta) i xi) / (1 - theta i xi)`,
  `xi = omega_h dt`; backward: its three-level recurrence).

It also fits the realized amplification factor, and evaluates the discrete
wave energy and the physical energy `0.5 sum V (beta p^2 + rho U^2)`.
Results: `postProcessing/cavityVerification/summary.json` and `history.csv`.
The backward discrete energy is reported only when every step is probed.

## Studies

```bash
./runStudy A     # scheme identification (6 runs)
./runStudy B     # temporal refinement (10 runs)
./runStudy C     # spatial refinement (5 runs)
./runStudy       # all three
```

Each run is a separate copy under `studies/`; `summarize_cavity_study.py`
writes the tables `studies/study*.csv` with observed orders.

`plot_cavity_results.py studies FIGURE.pdf TABLES.tex` writes the study-A
figure (modal amplitude and discrete-energy drift) and LaTeX table rows for
all three studies.
