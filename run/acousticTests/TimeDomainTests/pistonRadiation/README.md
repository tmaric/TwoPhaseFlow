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
mesh with `blockMesh`.
It then calls `setPMLFields`; `acousticWaveFoam` reads the generated `sigma`
field and derives its time-domain PML coefficients.

For a short smoke test:

```bash
MESH_CELLS_PER_LAMBDA=10 N_PERIODS=12 ./Allrun
```

The run reconstructs the latest pressure field and writes analytical
comparison data to `verification_axis.csv` and `verification_ring.csv`.
