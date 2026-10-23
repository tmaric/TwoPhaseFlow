# Time-domain piston radiation

Time-domain counterpart of `FrequencyDomainTests/pistonRadiation`, solved with
`acousticWaveFoam`. The cases share the same default drive, medium, piston,
domain, PML, and mesh parameters.

## Run

Load the OpenFOAM and TwoPhaseFlow environments, build the solver, then run:

```bash
wmake apps/utilities/mesh/setPMLFields
wmake solver/acousticWaveFoam
cd run/acousticTests/TimeDomainTests/pistonRadiation
./Allrun
```

Edit `caseParams.sh` to change the setup. `Allrun` calls `prepareCase`, which
renders the dependent OpenFOAM inputs before creating the rectangular wedge
mesh with `blockMesh`. The two-block mesh places an exact patch boundary at the
piston radius, so no subsequent face selection or patch reconstruction is
required.
It then calls `setPMLFields`; `acousticWaveFoam` reads the generated `sigma`
field and derives its time-domain PML coefficients.

For a short smoke test:

```bash
MESH_CELLS_PER_LAMBDA=10 N_PERIODS=12 ./Allrun
```

The run writes analytical comparison data directly from the parallel probe
histories to `verification_axis.csv` and `verification_ring.csv`. Automatic
field reconstruction is omitted because OpenFOAM's `reconstructPar` cannot
clone the time-dependent piston boundary condition reliably.

## Time integration and piston condition

`TIME_METHOD` selects the update: `blockCoupled` (default, monolithic);
`thetaNewmark` (segregated) and `fullyImplicit` (backward) are kept for
comparison only. The piston
boundary condition imposes `dp/dn = rho dU/dt` along the outward normal, the
same as the frequency case (`dP/dn = -i omega rho U`); before 2026-10-21 the
density factor was missing, which reduced the pressure by `1/rho`.

## Refinement study

`./runPistonStudy` runs the monolithic update at 20, 30 and 40 cells per
wavelength for 25 periods and writes `studies/studyPiston.csv`. Only the
on-axis comparison is meaningful: the ring comparison of
`verify_piston_radiation.py` uses the far-field (Fraunhofer) formula on a
ring inside the Rayleigh distance `a^2/lambda`, where it does not apply.
