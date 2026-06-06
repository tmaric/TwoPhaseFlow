# WedgeLevitatorDropAlpha case

This case is parameter-driven.  
Edit one file and regenerate all dependent inputs automatically.

## Files and roles

- `caseParams.sh`: single source of truth for drive, PML, mesh, solver, and drop setup.
- `prepareCase`: renders templates using values from `caseParams.sh`.
- `constant/transportProperties.in`: template for `constant/transportProperties`.
- `constant/levitatorWedgeHex.geo.in`: template for `constant/levitatorWedgeHex.geo`.
- `0.orig/Pim.in`: template for `0.orig/Pim`.
- `Allrun`: full workflow (clean, render, mesh, convert, set fields, MPI solve).
- `runHeightSweep.sh`: height sweep driver used when `RUN_MODE=fullFig4Sweep`.
- `plotHeightSweep.py`: plots the height-sweep CSV output.

## Typical usage

1. Edit parameters in:
   - `caseParams.sh`
2. Run full case:
   - `./Allrun`

`Allrun` already calls `./prepareCase`, so manual rendering is optional.

By default `RUN_MODE=singlePoint`, which preserves the original DropAlpha
single-case run. `Allrun` always initializes the ellipsoid with
`setAlphaField` and always solves with `acousticHelmholtzFoam`. To launch the
height sweep from this case:

```bash
RUN_MODE=fullFig4Sweep ./Allrun
```

For a small smoke test:

```bash
RUN_MODE=fullFig4Sweep SWEEP_POINTS=2 FIG4_ASPECT_RATIOS=1 RESULTS_DIR=sweepSmoke ./Allrun
```

## Optional: render-only step

Use this when you want to inspect generated files without running the solver:

```bash
./prepareCase
```

This updates:
- `constant/transportProperties`
- `constant/levitatorWedgeHex.geo`
- `0.orig/Pim`
- `system/setAlphaFieldDict`

## MPI run

`Allrun` runs decomposed MPI by default.  
Optional rank override:

```bash
NPROCS=8 ./Allrun
```

Drop initialization is controlled from `caseParams.sh`:
- `DROP_RADIUS`
- `DROP_HORIZONTAL_LONG_AXIS`
- `DROP_CENTER_Y`

Directional rectangular-PML damping is controlled with:
- `SIGMA_MAX_X`
- `SIGMA_MAX_Y`
- `SIGMA_MAX_Z`

`prepareCase` writes these values to `sigmaMaxVector`. The scalar
`SIGMA_MAX` remains available for backward compatibility and initializes the
x- and y-direction values when they are not set explicitly.

## Environment

Run with your OpenFOAM environment loaded:

```bash
source ~/TwoPhaseFlow/scripts/bashrc
source ~/openfoam/etc/bashrc
```

Then execute:

```bash
cd /home/minkowski/TwoPhaseFlow/run/acousticTests/FrequencyDomainTests/Andrade2019/WedgeLevitatorDropAlpha
./Allrun
```
