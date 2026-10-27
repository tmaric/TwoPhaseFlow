# layeredInterface1D

Heterogeneous/interface validation case for `acousticHelmholtzFoam`.

## What this case checks
- Normal-incidence acoustic reflection/transmission at a planar material interface.
- Material jump is represented with `alpha.water`:
  - `alpha.water = 0`: medium 1 (`rhog`, `cg`)
  - `alpha.water = 1`: medium 2 (`rhol`, `cl`)
- The analytical pressure covers the gas, liquid, and right-end PML.
- The reported complex-pressure error is evaluated over the complete domain.

## Run (MPI default)
```bash
./Allrun
```

## Outputs
- Solver log: `log.acousticHelmholtzFoam`
- Verification report: `postProcessing/interfaceValidation/1/metrics.txt`
- Component comparison: `postProcessing/interfaceValidation/1/pressureField_Pre_Pim_compare.png`
- Amplitude comparison: `postProcessing/interfaceValidation/1/pressureField_abs_compare.png`
- Comparison data: `postProcessing/interfaceValidation/1/pressureFieldComparison.csv`

## PML damping sensitivity

The damping sweep holds the mesh, PML thickness, and polynomial order fixed:

```bash
./run_pml_sensitivity.sh
```

Outputs are written to `postProcessing/pmlSensitivity`.

## Mesh sensitivity

The mesh sweep keeps the material interface and PML start aligned with mesh
faces and uses the same whole-domain complex-pressure error:

```bash
./run_pml_mesh_sensitivity.sh
```

Outputs are written to `postProcessing/pmlMeshSensitivity`.

## Face-area comparison

Both material orderings at `NX = 8000`, the five damping values and the five
aligned meshes are solved with `areaFraction legacy` (original donor-based
selection) and `areaFraction plicAverage` (geometric averaging), both with
`flux legacy`, serially and on 4 and 8 ranks:

```bash
python3 run_face_area_comparison.py --output studyResults/faceAreaComparison
```

`--configs` (`NX:sigmaMax:forward|reverse`) and `--ranks` select a subset.
`face-area-comparison.csv` reports `referenceDifference`, the relative
difference from the serial `legacy` pressure of the same configuration,
`parallelDifference` from the serial pressure of the same mode, and
`pressureRelL2` at the cell centres. The sweeps report `P_relL2` at 1200
points interpolated from the cell values, so the two errors differ.

## Interpretation
`metrics.txt` reports `P_relL2`, the relative complex-pressure error over all
uniformly sampled points from the driven boundary through the outer PML
boundary. Component and pressure-amplitude errors are retained as secondary
diagnostics; the paper uses `P_relL2` as the quantitative measure.
