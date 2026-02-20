# pistonRadiation case

This case is parameter-driven.  
Edit one file and regenerate all dependent inputs automatically.

## Files and roles

- `caseParams.sh`: single source of truth for frequency, PML, and mesh controls.
- `prepareCase`: renders templates using values from `caseParams.sh`.
- `constant/transportProperties.in`: template for `constant/transportProperties`.
- `constant/pistonRadiation.geo.in`: template for `constant/pistonRadiation.geo`.
- `0.orig/Pim.in`: template for `0.orig/Pim`.
- `Allrun`: full workflow (clean, render, mesh, convert, patch, solve).

## Typical usage

1. Edit parameters:
   - `caseParams.sh`
2. Run full case:
   - `./Allrun`

`Allrun` already calls `./prepareCase`, so manual rendering is not required.

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
