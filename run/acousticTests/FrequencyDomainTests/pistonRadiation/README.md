# pistonRadiation case

This case is parameter-driven.  
Edit one file and regenerate all dependent inputs automatically.

## Files and roles

- `caseParams.sh`: single source of truth for frequency, PML, and mesh controls.
- `prepareCase`: renders templates using values from `caseParams.sh`.
- `constant/transportProperties.in`: template for `constant/transportProperties`.
- `constant/pistonRadiation.geo.in`: template for `constant/pistonRadiation.geo`.
- `templates/Pim.in`: template for `0.orig/Pim`.
- `Allrun`: switchable workflow for serial reference or MPI solver.

## Typical usage

1. Edit parameters:
   - `caseParams.sh`
2. Run full case:
   - `./Allrun`

`Allrun` already calls `./prepareCase`, so manual rendering is not required.

## Solver execution

`Allrun` supports both modes:

```bash
# serial reference (default)
./Allrun
# or explicitly
RUN_MODE=serial ./Allrun

# MPI/decomposed run
RUN_MODE=mpi ./Allrun
# optional override for MPI ranks
NPROCS=8 RUN_MODE=mpi ./Allrun
```

Default mode is set in `caseParams.sh` via `RUN_MODE`.

Serial log file:
- `log.freBlockCoupledFoam`

MPI log file:
- `log.freBCMFoam`

## Optional: render-only step

Use this when you want to inspect generated files without running the solver:

```bash
./prepareCase
```

This updates:
- `constant/transportProperties`
- `constant/pistonRadiation.geo`
- `0.orig/Pim`

## Environment

Run with your usual OpenFOAM environment loaded, for example:

```bash
source ~/TwoPhaseFlow/scripts/bashrc
source ~/openfoam/etc/bashrc
```

Then execute:

```bash
# from this case directory
./Allrun
```

## Analytical comparison (COMSOL-style on-axis)

After a solve, run:

```bash
python3 postprocess_compare.py
```

Optional far-field fit radii (inside non-PML region):

```bash
python3 postprocess_compare.py --r-src 0.2 --r-far 4.0
```

Outputs are written to:

- `postProcessing/analyticalCompare/<latestTime>/onAxisComparison.png`
- `postProcessing/analyticalCompare/<latestTime>/onAxisComparison.csv`
- `postProcessing/analyticalCompare/<latestTime>/farFieldPatternComparison.png`
- `postProcessing/analyticalCompare/<latestTime>/farFieldPatternComparison.csv`
- `postProcessing/analyticalCompare/<latestTime>/metrics.txt`
