# WedgeLevitatorHeightSweep

Frequency-domain empty-levitator sweep with serial/MPI switch.

This case reuses the parameterized setup of `WedgeLevitatorDropAlpha`, but with
an empty cavity (`alpha.water = 0`) and a sweep in reflector height ratio:

- `D/lambda` from `0.5` to `1.5`
- `100` sample points

## Files

- `caseParams.sh`: base parameters, run mode (`RUN_MODE`), and sweep controls.
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

Mode switch:

```bash
RUN_MODE=serial ./Allrun
RUN_MODE=mpi ./Allrun
```

## Run full sweep + plot

```bash
./runHeightSweep.sh
```

MPI sweep:

```bash
RUN_MODE=mpi ./runHeightSweep.sh
```

Outputs:

- `sweepResults/heightSweepResults.csv`
- `sweepResults/heightSweep_force.png`

CSV columns:

- `index`: sweep index `[0..SWEEP_POINTS-1]`
- `height_mm`: geometric gap in millimeters (plot x-axis)
- `D_m`: geometric gap in meters
- `height_fac`: `D/lambda`
- `Fn_reflector1_wedge_N`: force from `pr` integration on patch `reflector1` on the wedge sector
- `Fn_reflector1_axisym_N`: axisymmetric full-360deg force (`Fn_reflector1_wedge_N * 360/WEDGE_DEG`)
- `axisym_factor`: `360/WEDGE_DEG`
