# WedgeLevitatorDynamicDrop case

This case is parameter-driven.  
Edit one file and regenerate all dependent inputs automatically.

## Files and roles

- `caseParams.sh`: single source of truth for drive, PML, and mesh controls.
  - includes `DROP_RADIUS` and `DROP_CENTER_Y` for initial drop size/vertical position (`setAlphaField`).
- `prepareCase`: renders templates using values from `caseParams.sh`.
- `constant/transportProperties.in`: template for `constant/transportProperties`.
- `constant/levitatorWedgeHex.geo.in`: template for `constant/levitatorWedgeHex.geo`.
- `0.orig/Pim.in`: template for `0.orig/Pim`.
- `Allrun`: full workflow with decomposed MPI run.

## Typical usage

1. Load environment (required):
   - `source ~/TwoPhaseFlow/scripts/bashrc`
   - `source ~/openfoam/etc/bashrc`
2. Edit parameters in:
   - `caseParams.sh`
3. Run full case:
   - `./Allrun`

`Allrun` already calls `./prepareCase`, so manual rendering is optional.
`Allrun` also checks `gmsh` and auto-builds the selected solver if needed.

## Contributor quick check (recommended)

Run these once before first case run:

```bash
source ~/TwoPhaseFlow/scripts/bashrc
source ~/openfoam/etc/bashrc
cd /home/minkowski/TwoPhaseFlow
wmake solver/interFALFlow
```

Then in this case directory:

```bash
./Allrun
```

## MPI run

`Allrun` always runs decomposed MPI:

```bash
./Allrun
# optional override for MPI ranks
NPROCS=8 ./Allrun
```

Log file:
- `log.interFALFlow`

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

## Environment

Run with your OpenFOAM environment loaded:

```bash
source ~/TwoPhaseFlow/scripts/bashrc
source ~/openfoam/etc/bashrc
```

Then execute:

```bash
cd /home/minkowski/TwoPhaseFlow/run/acousticTests/FrequencyDomainTests/Andrade2019/WedgeLevitatorDynamicDrop
./Allrun
```
