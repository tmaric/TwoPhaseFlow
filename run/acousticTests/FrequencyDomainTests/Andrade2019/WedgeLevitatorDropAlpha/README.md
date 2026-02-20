# WedgeLevitatorDropAlpha case

This case is parameter-driven.  
Edit one file and regenerate all dependent inputs automatically.

## Files and roles

- `caseParams.sh`: single source of truth for drive, PML, and mesh controls.
- `prepareCase`: renders templates using values from `caseParams.sh`.
- `constant/transportProperties.in`: template for `constant/transportProperties`.
- `constant/levitatorWedgeHex.geo.in`: template for `constant/levitatorWedgeHex.geo`.
- `0.orig/Pim.in`: template for `0.orig/Pim`.
- `Allrun`: full workflow (clean, render, mesh, convert, set fields, solve).

## Typical usage

1. Edit parameters in:
   - `caseParams.sh`
2. Run full case:
   - `./Allrun`

`Allrun` already calls `./prepareCase`, so manual rendering is optional.

## Optional: render-only step

Use this when you want to inspect generated files without running the solver:

```bash
./prepareCase
```

This updates:
- `constant/transportProperties`
- `constant/levitatorWedgeHex.geo`
- `0.orig/Pim`

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
