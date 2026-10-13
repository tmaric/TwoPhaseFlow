# twoLayerCavity1D

Free eigenmode of a closed rigid one-dimensional cavity with a gas layer
on `[0, x_I]` and a liquid layer on `[x_I, L]`, solved with
`acousticWaveFoam`. There is no PML, no boundary forcing and no added
damping, so the case tests the interface discretization, the velocity
recovery beside the interface, and parallel consistency.

## Configuration

- Air (`RHOG`, `CG`) and water (`RHOL`, `CL`); rigid walls at both ends.
- Layer lengths in units `u = L_GAS/GAS_UNITS`. The defaults
  (`L_GAS = 0.343 m`, 32 and 35 units) place the second gas half-wave
  resonance and the first liquid quarter-wave resonance near 1 kHz, so the
  two modes there hybridize and store energy in both layers
  (994.517 Hz and 1005.067 Hz).
- `INTERFACE_OFFSET = s` places the interface on a face (`s = 0`) or inside
  a cell (`0 < s < 1`), at the same relative position on every mesh.
- `MODE = lower | upper` selects the exact mode nearest to `F_TARGET`.
  The initial pressure is that mode at the cell centres, with
  `pDot = U = 0`; `alpha.water` is the exact volume fraction.
- `NPROCS = 2` decomposes manually at the interface: the pure gas cells form
  rank 0, so the processor face is the interface face, or the face beside
  the interface cell.

## Run

Load OpenFOAM-v2606 and `scripts/bashrc`, build `acousticWaveFoam`, then

```bash
./Allrun
```

All entries of `caseParams.sh` can be overridden from the environment, e.g.

```bash
MODE=upper INTERFACE_OFFSET=0.25 CELLS_PER_UNIT=4 ./Allrun
```

`prepareCase` calls `prepare_two_layer.py`, which solves the exact
frequency equation (`twoLayerModes.py`) and writes all run-dependent inputs.

## Post-processing

`postprocess_two_layer.py` builds the discrete reference from the face
coefficients the solver actually used (`invRhof`, written by `writeObjects`)
and the cell coefficients of `alpha.water`. It projects the probed pressure
onto the discrete mode nearest to the exact frequency with the left
eigenvector, so other modes in the initial field do not affect the
projected amplitude. Unequal coefficients on the two sides of a processor
face are reported and make the operator nonsymmetric.

It reports the deviation from the fully discrete recurrence, the fitted
frequency and decay rate, the frequency error of the discrete mode, the
pressure and velocity errors against the exact mode, and the discrete
energy. Results: `postProcessing/twoLayerVerification/`.

## Studies

```bash
./runStudy A     # scheme identification (8 runs)
./runStudy B     # spatial refinement, interface on a face (10 runs)
./runStudy C     # spatial refinement, interface within a cell (20 runs)
./runStudy D     # serial versus two ranks (8 runs)
./runStudy       # all, about 12 minutes on one core
```

`summarize_two_layer_study.py` writes `studies/study*.csv`;
`plot_two_layer_results.py studies FIGURE.pdf TABLES.tex` writes the
frequency-convergence figure and LaTeX table rows.

## Radiation averages

`pr` and `momFlux` are probed with the fields. With `RADIATION_CHECK=1`
(serial, interface on a face), `postprocess_radiation.py` compares the
running averages over `RADIATION_PERIODS` periods with the fully discrete
solution, expanded in all discrete modes of the assembled operator, the
exact mode averaged over the same window, and the exact means, separately
for the gas and liquid interiors and the two cells next to the interface.
`./runRadiationStudy` refines both hybrid modes (8 runs) and writes
`studies/radiation/studyC.csv`.
