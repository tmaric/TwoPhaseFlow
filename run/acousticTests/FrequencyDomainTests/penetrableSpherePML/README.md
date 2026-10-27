# penetrableSpherePML case

Curved transmitting interface with PML truncation. A fluid sphere of radius
`a = 1 mm` lies on the axis of a one-degree wedge and is represented only by
`setAlphaField` volume fractions. An axial point source at `y = 2.5a` is an
excluded box whose faces carry the exact total pressure. The physical region
`r <= 4a, -4a <= y <= 6a` has uniform cells `h = a/N` and exterior `ka = 1`.

## Files and roles

- `sphere_pml.py`: generates, runs and evaluates the cases (`blockMesh`,
  `topoSet`, `subsetMesh`, `setAlphaField`, `setPMLFields`,
  `acousticHelmholtzFoam`). Uses the binaries in `PATH`.
- `sphere_reference.py`: exact pressure of the point source and the fluid
  sphere (multipole series, 61 orders). `python3 sphere_reference.py` runs its
  checks: interface conditions with analytic derivatives, the incident
  expansion, the homogeneous limit and the series truncation.
- `summarize_sphere_pml.py`: writes `convergence.csv` per study and the
  serial versus eight-rank comparison.
- `Allrun`: the full study; `Allclean` removes `studyResults/`.

The same sphere with exact boundary data and no PML is `../penetrableSphere`.

## Truncations

- `dirichlet`: exact total pressure on the outer boundary of the physical
  region, no PML.
- `pml`: Cartesian PML of thickness `L` around the physical region with
  `wedgeAzimuthalStretch true`, i.e. the cylindrical radial PML.
- `pmlCartesian`: the same PML with `wedgeAzimuthalStretch false`.

The PML strength is given as `sigmaMax/omega` (`--sigma-factor`) and its
thickness in units of `a` (`--pml-thickness`). The interface method is
`areaFraction plicAverage; flux legacy;`.

## Material pairs

| Name | Host rho, c | Drop rho, c |
|---|---|---|
| `waterInAir` | 1.2, 343 | 1000, 1500 |
| `oilInWater` | 1000, 1500 | 920, 1450 |
| `airInWater` | 1000, 1500 | 1.2, 343 |

The host occupies the gas slot (`rhog`, `cg`) of the solver.

## Usage

Load the OpenFOAM and TwoPhaseFlow environment, then

```bash
./Allrun
```

`SPHERE_SIZES` (default `8,16,32,64`) and `SPHERE_RANK_SIZES` (default
`16,32`) select the resolutions. A single study:

```bash
python3 sphere_pml.py --output studyResults/smoke \
    --materials waterInAir --sizes 8,16
```

## Output

Each study directory under `studyResults/` contains `settings.json`
(arguments, Git revision, binary and script hashes), `results.json`, and per
material `reference_check.json`, the prepared cases `prepared-<truncation>-N<N>`
and the solved cases `<truncation>-N<N>`. Errors are volume-weighted relative
complex-pressure L2 norms in the physical region: total, host, drop, interface
cells and far host field (`r > 2.5a`, outside the source neighbourhood).
`minusDirichletRelL2` is the difference between a PML run and the `dirichlet`
run on the same physical mesh.
