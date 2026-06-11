# WedgeLevitatorDropAlpha case

This case is parameter-driven.  
Edit one file and regenerate all dependent inputs automatically.

## Files and roles

- `caseParams.sh`: single source of truth for drive, PML, mesh, solver, and drop setup.
- `prepareCase`: renders templates using values from `caseParams.sh`.
- `constant/transportProperties.in`: template for `constant/transportProperties`.
- `generateCfMeshInputs.py`: generates the cfMesh surface and dictionaries.
- `0.orig/Pim.in`: template for `0.orig/Pim`.
- `Allrun`: full workflow (clean, render, cfMesh, set fields, MPI solve).
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
- `constant/levitator2D.fms`
- `0.orig/Pim`
- `system/setAlphaFieldDict`
- `system/meshDict`
- `system/extrudeMeshDict`

The case uses cfMesh. After `setAlphaField`, the workflow selects
mixed cells satisfying `0 < alpha.water < 1`, optionally grows the set by
face-adjacent layers on both liquid and gas sides, refines the selected cells,
and reconstructs `alpha.water`. Configure it in `caseParams.sh` using:

- `CFMESH_MAX_CELL_SIZE`
- `ALPHA_INTERFACE_REFINEMENT`
- `ALPHA_INTERFACE_REFINEMENT_LEVELS`
- `ALPHA_INTERFACE_ADJACENT_LAYERS`

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

## Performance studies

`runPerformanceScaling.sh` generates CSV data for strong scaling and
problem-size scaling. Strong scaling prepares one mesh and reuses it across
MPI-rank counts. Problem-size scaling varies `CFMESH_MAX_CELL_SIZE` while
holding the MPI-rank count fixed.

```bash
# Run both studies with the defaults
./runPerformanceScaling.sh

# Fixed mesh with selected MPI-rank counts
STRONG_RANKS="1 2 4 8 16" ./runPerformanceScaling.sh strong

# Selected mesh sizes at eight MPI ranks
SIZE_CELL_SIZES="8e-5 6e-5 4e-5 3e-5" SIZE_RANKS=8 \
    ./runPerformanceScaling.sh size
```

Results and per-run logs are written under `performanceResults/`. The CSV files
report the cell count, block-system unknown count (`2*cells`), measured elapsed
time, and solver-reported clock time. The strong-scaling CSV additionally
reports speedup and parallel efficiency relative to the first rank count. A
single unreported warm-up solve prevents one-time coded-function compilation
from contaminating the timings. Set `PERFORMANCE_WARMUP=0` to disable it.
Additional launcher arguments can be supplied with `MPI_OPTIONS`; OpenMPI uses
`--oversubscribe` by default to match the standard OpenFOAM run helper.

Generate publication-ready PNG and PDF plots with:

```bash
# Plot the newest benchmark under performanceResults/
./plotPerformanceScaling.py

# Plot a selected benchmark directory
./plotPerformanceScaling.py performanceResults/20260611-120000
```

The figures are written to `<results_dir>/figures/`. Strong-scaling output
contains runtime, speedup, and efficiency panels. Problem-size output contains
total runtime and runtime per block-system unknown.

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
