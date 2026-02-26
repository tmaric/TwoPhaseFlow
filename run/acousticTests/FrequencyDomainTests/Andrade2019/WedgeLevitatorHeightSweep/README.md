# WedgeLevitatorHeightSweep

Frequency-domain empty-levitator sweep for `freBlockCoupledFoam`.

This case reuses the parameterized setup of `WedgeLevitatorDropAlpha`, but with
an empty cavity (`alpha.water = 0`) and a sweep in reflector height ratio:

- `D/lambda` from `0.5` to `1.5`
- `100` sample points

## Files

- `caseParams.sh`: base parameters and sweep controls.
- `prepareCase`: renders `constant/transportProperties`, `constant/levitatorWedgeHex.geo`, and `0.orig/Pim`.
- `Allrun`: runs one single height point (value from `HEIGHT_FAC`).
- `runHeightSweep.sh`: executes all 100 height points and collects force data.
- `plotHeightSweep.py`: makes the plot from CSV.

## Run one point

```bash
# default HEIGHT_FAC from caseParams.sh
./Allrun

# override one point, e.g. D/lambda = 0.8
HEIGHT_FAC=0.8 ./Allrun
```

## Run full sweep + plot

```bash
./runHeightSweep.sh
```

Outputs:

- `sweepResults/heightSweepResults.csv`
- `sweepResults/heightSweep_force.png`

CSV columns:

- `index`: sweep index `[0..SWEEP_POINTS-1]`
- `height_fac`: `D/lambda`
- `D_m`: geometric gap in meters
- `Fn_reflector1_N`: normal force from `pr` integration on patch `reflector1`
