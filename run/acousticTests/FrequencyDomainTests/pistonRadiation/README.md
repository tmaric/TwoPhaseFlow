# pistonRadiation case

This case is parameter-driven.  
Edit one file and regenerate all dependent inputs automatically.

## Files and roles

- `caseParams.sh`: single source of truth for frequency, PML, and mesh controls.
- `prepareCase`: renders templates using values from `caseParams.sh`.
- `constant/transportProperties.in`: template for `constant/transportProperties`.
- `system/blockMeshDict.in`: template for the rectangular wedge mesh.
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
The mesh is generated directly with `blockMesh`.
The piston edge is an exact radial block boundary at `PISTON_RADIUS`; it does
not move when the cells-per-wavelength setting changes.

The meridional domain is rectangular.  The physical region extends to
`PML_RMIN` in the radial and axial directions, and the PML continues to
`PML_RMAX`.  Damping is active only towards the outer radial and upper
boundaries; the symmetry axis and rigid baffle are not damped.

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
- `system/blockMeshDict`
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

The default exterior-field reconstruction follows the COMSOL baffled-piston
setup: a spherical source surface at PML_RMIN is completed across the rigid
baffle using sound-hard symmetry and evaluated with the full finite-distance
Kirchhoff--Helmholtz integral.

Optional reconstruction settings include changing the source radius, disabling
the symmetry completion for diagnostics, and running the analytical quadrature
validation:

```bash
python3 postprocess_compare.py --r-src 0.18 --r-far 4.0
python3 postprocess_compare.py --source-symmetry none --output-name hemisphereDiagnostic
python3 validate_kirchhoff.py
```

Outputs are written to:

- `postProcessing/analyticalCompare/<latestTime>/onAxisComparison.png`
- `postProcessing/analyticalCompare/<latestTime>/onAxisComparison.csv`
- `postProcessing/analyticalCompare/<latestTime>/farFieldPatternComparison.png`
- `postProcessing/analyticalCompare/<latestTime>/farFieldPatternComparison.csv`
- `postProcessing/analyticalCompare/<latestTime>/metrics.txt`

## Mesh convergence batch

Run the near- and far-field analytical comparisons for 20, 30, 40, 60, and
80 cells per wavelength:

```bash
./meshConv.sh
```

Outputs are archived in `meshConvergence/`, including per-resolution comparison
CSVs/plots, `metrics_summary.csv`, aggregate convergence plots, and LaTeX-ready
tables:

- `meshConvergence/nearField_convergence_table.tex`
- `meshConvergence/farField_convergence_table.tex`

The summary CSV reports `h_over_lambda`, near-field relative `L2` and `Linf`
errors, and the far-field pressure-magnitude relative `L2` error.

Optional overrides:

```bash
MESH_CONV_RESOLUTIONS="20 50" ./meshConv.sh
POSTPROCESS_ARGS="--n-points 300 --n-theta 91" ./meshConv.sh
```
