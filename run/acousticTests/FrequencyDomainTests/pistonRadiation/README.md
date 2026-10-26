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
`wedgeAzimuthalStretch true` in `constant/transportProperties` makes
`setPMLFields` add the azimuthal damping `sigma.zz = (1/r) int sigma.xx dr`,
so that the Cartesian PML equals the radial (cylindrical) PML in the wedge.

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

## Comparison with the exact field

After a solve, run:

```bash
python3 postprocess_compare.py
```

The reference is the exact pressure of a uniform piston in a rigid baffle (the
Rayleigh integral), evaluated in `piston_reference.py` as plane wave minus edge
wave: the integral over the piston surface reduces exactly to a line integral
over the rim.  On the axis this gives the closed form
`rho c u0 [exp(ikz) - exp(ik sqrt(z^2 + a^2))]`.  `validate_reference.py` checks
the reference against the on-axis closed form, a direct surface quadrature of
the Rayleigh integral, node doubling and continuity at the rim:

```bash
python3 validate_reference.py
```

The cell-centre pressure is compared with the reference at the same cell
centres, without interpolation and without an extrapolation surface:

- `field_relL2`: volume-weighted complex-pressure error in the physical region
  (`r, z <= PML_RMIN`);
- `axis_relL2`, `axis_relLinf`: complex-pressure errors in the cell row next to
  the symmetry axis;
- `*_magnitude_relL2`: the same norms of `|p|`.

The pressure level on the arc `r = 0.9 PML_RMIN` (`--r-arc`) is probed for
plotting only.  Outputs are written to
`postProcessing/analyticalCompare/<latestTime>/`:

- `onAxisComparison.csv`, `onAxisComparison.png`
- `arcComparison.csv`, `arcComparison.png`
- `fieldError.png` (pointwise error in the physical region)
- `fieldComparison.csv.gz` (all cells of the physical region)
- `metrics.txt`

## Mesh convergence batch

Run the comparison for 20, 30, 40, 60, and 80 cells per wavelength:

```bash
./meshConv.sh
```

`meshConv.sh` first runs `validate_reference.py`.  Outputs are archived in
`meshConvergence/`: per-resolution comparison CSVs and plots,
`metrics_summary.csv`, aggregate plots and the LaTeX table
`meshConvergence/piston_convergence_table.tex` (errors and observed orders).

Optional overrides:

```bash
MESH_CONV_RESOLUTIONS="20 50" ./meshConv.sh
POSTPROCESS_ARGS="--n-theta 91" ./meshConv.sh
```
