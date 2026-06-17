# pistonRadiation case

This case is parameter-driven.  
Edit one file and regenerate all dependent inputs automatically.

## Files and roles

- `caseParams.sh`: single source of truth for frequency, PML, and mesh controls.
- `prepareCase`: renders templates using values from `caseParams.sh`.
- `constant/transportProperties.in`: template for `constant/transportProperties`.
- `constant/pistonRadiation.geo.in`: template for `constant/pistonRadiation.geo`.
- `templates/Pim.in`: template for `0.orig/Pim`.
- `Allrun`: workflow for decomposed MPI solver.

## Typical usage

1. Load environment (required):
   - `source ~/TwoPhaseFlow/scripts/bashrc`
   - `source ~/openfoam/etc/bashrc`
2. Edit parameters:
   - `caseParams.sh`
3. Run full case:
   - `./Allrun`

`Allrun` already calls `./prepareCase`, so manual rendering is not required.
`Allrun` also checks `gmsh` and auto-builds the solver if the executable is missing.

## Contributor quick check (recommended)

Run these once before first case run:

```bash
source ~/TwoPhaseFlow/scripts/bashrc
source ~/openfoam/etc/bashrc
cd /home/minkowski/TwoPhaseFlow
wmake solver/acousticHelmholtzFoam
```

Then in this case directory:

```bash
./Allrun
```

## Solver execution

`Allrun` runs decomposed MPI by default:

```bash
./Allrun
```

Set the MPI rank count with `numberOfSubdomains` in
`system/decomposeParDict`.

Log file:
- `log.acousticHelmholtzFoam`

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

## Mesh convergence batch

Run near-field/on-axis and far-field analytical comparisons for 20, 50, and
100 cells per wavelength:

```bash
./meshConv.sh
```

Outputs are archived in `meshConvergence/`, including per-resolution comparison
CSVs/plots, `metrics_summary.csv`, aggregate convergence plots, and LaTeX-ready
tables:

- `meshConvergence/nearField_convergence_table.tex`
- `meshConvergence/farField_convergence_table.tex`

The summary CSV reports `h_over_lambda`, near-field relative `L2` and `Linf`
errors, far-field relative `L2` and `Linf` errors, and the observed convergence
orders between successive mesh resolutions.

Optional overrides:

```bash
MESH_CONV_RESOLUTIONS="20 50" ./meshConv.sh
POSTPROCESS_ARGS="--n-points 300 --n-theta 91" ./meshConv.sh
```
