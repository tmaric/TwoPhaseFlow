# WedgeLevitatorHeightSweep

Frequency-domain height sweep for reproducing Fig. 4 of Andrade et al. (2019).

The case sweeps reflector height and can initialize either:

- an empty levitator, or
- a constant-volume spheroid with equivalent spherical radius `1 mm`,
  centered at `H/2`.

The default full sweep reproduces the simulation curves from Fig. 4:

- empty levitator
- sphere (`a/b = 1`)
- oblate spheroid (`a/b = 2`)
- oblate spheroid (`a/b = 3`)

## Files

- `caseParams.sh`: base parameters, solver, and sweep controls.
- `prepareCase`: renders `constant/transportProperties`, `constant/levitatorWedgeHex.geo`, `0.orig/Pim`, and `system/setExprFieldsDict`.
- `Allrun`: runs one single height point (value from `HEIGHT_FAC`) with either empty or spheroid initialization.
- `runHeightSweep.sh`: executes the Fig. 4 shape set over all height points and collects force data.
- `plotHeightSweep.py`: makes the plot from CSV.

## Run directly

With the default `RUN_MODE=fullFig4Sweep` in [caseParams.sh](/home/local/CSI/cx80jevu/TwoPhaseFlow/run/acousticTests/FrequencyDomainTests/Andrade2019/WedgeLevitatorHeightSweep/caseParams.sh), a plain

```bash
./Allrun
```

launches the full Fig. 4 sweep and writes the reproduced figure.

If you want `./Allrun` to run only a single case, switch `RUN_MODE=singlePoint`.

## Run one point

```bash
# after setting RUN_MODE=singlePoint in caseParams.sh
./Allrun

# override one point, e.g. D/lambda = 0.8 with a/b = 3
DROP_MODE=oblate DROP_ASPECT_RATIO=3 HEIGHT_FAC=0.8 ./Allrun

# sphere at the cavity mid-plane
DROP_MODE=oblate DROP_ASPECT_RATIO=1 ./Allrun

# override one point, empty levitator
HEIGHT_FAC=0.8 ./Allrun
```

MPI run:

```bash
NPROCS=8 ./Allrun
```

## Run full Fig. 4 sweep + plot

```bash
./runHeightSweep.sh
```

Optional rank override for sweep:

```bash
NPROCS=8 ./runHeightSweep.sh
```

Outputs:

- `sweepResults/heightSweepResults.csv`
- `sweepResults/reproducedFig4.png`
- `sweepResults/reproducedFig4.pdf`

The output directory and figure basename are controlled by
`RESULTS_DIR` and `FIG4_OUTPUT_BASENAME` in
[caseParams.sh](/home/local/CSI/cx80jevu/TwoPhaseFlow/run/acousticTests/FrequencyDomainTests/Andrade2019/WedgeLevitatorHeightSweep/caseParams.sh).

CSV columns:

- `series_key`: curve identifier used by the plotter
- `series_label`: legend label
- `drop_mode`: `empty` or `oblate`
- `aspect_ratio`: `a/b`
- `index`: sweep index `[0..SWEEP_POINTS-1]`
- `height_mm`: geometric gap in millimeters (plot x-axis)
- `D_m`: geometric gap in meters
- `height_fac`: `D/lambda`
- `Fn_reflector1_wedge_N`: force from `pr` integration on patch `reflector1` on the wedge sector
- `Fn_reflector1_axisym_N`: axisymmetric full-360deg force (`Fn_reflector1_wedge_N * 360/WEDGE_DEG`)
- `axisym_factor`: `360/WEDGE_DEG`
