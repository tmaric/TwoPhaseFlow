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

## Interpretation
`metrics.txt` reports `P_relL2`, the relative complex-pressure error over all
uniformly sampled points from the driven boundary through the outer PML
boundary. Component and pressure-amplitude errors are retained as secondary
diagnostics; the paper uses `P_relL2` as the quantitative measure.
