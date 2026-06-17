# WedgeLevitatorDropContactDynamic case

This case is parameter-driven.
Edit one file and regenerate all dependent inputs automatically.

This variant starts from a 1 mm water drop attached to the reflector.
The initial sphere center is placed on the reflector plane (`DROP_CENTER_Y=0`),
which gives a hemispherical initial drop, and `reflector1` uses a
`dynamicAlphaContactAngle` boundary condition with equilibrium angle
`theta0 90`, advancing angle `thetaA 110`, receding angle `thetaR 70`,
and velocity scale `uTheta 0.01`.

## Files and roles

- `caseParams.sh`: single source of truth for drive, PML, and mesh controls.
  - includes `DROP_RADIUS` and `DROP_CENTER_Y` for initial drop size/vertical position (`setAlphaField`).
- `prepareCase`: renders templates using values from `caseParams.sh`.
  - generates `constant/levitator2D.fms`, `system/meshDict`, and `system/extrudeMeshDict`.
- `constant/transportProperties.in`: template for `constant/transportProperties`.
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
`Allrun` uses cfMesh and auto-builds the selected solver if needed.

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
```

Set the MPI rank count with `numberOfSubdomains` in
`system/decomposeParDict`.

Log file:
- `log.interFALFlow`
- Drop geometry tracking:
  - `postProcessing/dropGeometry/geometry.dat`
  - columns: `time`, `centerY`, `horizontalAxis`, `verticalAxis`, `xMax`, `yMin`, `yMax`, `weightedVolume`, `selectedCells`
  - plot script:
    - `python3 plotDropGeometry.py`
    - outputs `postProcessing/dropGeometry/dropGeometry.png` and `.pdf`

## Optional: render-only step

Use this when you want to inspect generated files without running the solver:

```bash
./prepareCase
```

This updates:
- `constant/transportProperties`
- `0.orig/Pim`
- `system/setAlphaFieldDict`
- `constant/levitator2D.fms`
- `system/meshDict`
- `system/extrudeMeshDict`

## Environment

Run with your OpenFOAM environment loaded:

```bash
source ~/TwoPhaseFlow/scripts/bashrc
source ~/openfoam/etc/bashrc
```

Then execute:

```bash
cd /home/minkowski/TwoPhaseFlow/run/acousticTests/FrequencyDomainTests/Andrade2019/WedgeLevitatorDropContactDynamic
./Allrun
```
