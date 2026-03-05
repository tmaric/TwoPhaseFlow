# layeredInterface1D

Heterogeneous/interface validation case for `acousticHelmholtzFoam`.

## What this case checks
- Normal-incidence acoustic reflection/transmission at a planar material interface.
- Material jump is represented with `alpha.water`:
  - `alpha.water = 0`: medium 1 (`rhog`, `cg`)
  - `alpha.water = 1`: medium 2 (`rhol`, `cl`)
- Right-end PML damps transmitted wave to approximate semi-infinite medium 2.

## Run
```bash
./Allrun
```

## Outputs
- Solver log: `log.acousticHelmholtzFoam`
- Verification report: `postProcessing/interfaceValidation/1/metrics.txt`
- Plot: `postProcessing/interfaceValidation/1/interface_RT_compare.png`

## Interpretation
`metrics.txt` reports numerical vs analytical complex reflection/transmission coefficients:
- `R = B/A`
- `T = C/A`

Analytical references (normal incidence):
- `R_th = (Z2 - Z1)/(Z2 + Z1)`
- `T_th = 2*Z2/(Z2 + Z1)`
with impedances `Zi = rhoi*ci`.
