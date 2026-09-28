# Sphere-normal reconstruction ladder

Static convergence of the PLIC interface normal reconstructed by `plicRDF`,
`gradAlpha` and `isoAlpha` from the volume fraction of a sphere (R = 1 mm in a
4R box, centre off every mesh symmetry plane) on uniform hex meshes with
N = 20, 40, 80, 160 cells per side (R/h = 5, 10, 20, 40). One reconstruction
step per case (`reconstructInterface`, `stopAt writeNow`).

The error columns were appended to `reconstructionError` (`src/postProcessing`)
for this benchmark, after the existing twelve, which are unchanged:

| column | meaning |
|---|---|
| `LNormalDiff1`, `LNormalDiff2`, `LNormalDiffInf` | mean, rms, max over the interface cells (`|normal| != 0`) of `|n - n_exact|`, `n` the unit PLIC normal, `n_exact` the sphere normal at the reconstructed centre (both point into alpha = 1) |
| `areaSum` | `sum |normal|` over the cells: `normal_` is the PLIC area vector, so this is the reconstructed interface area |
| `areaExact` | `4 pi R^2` for `type sphere`, 0 otherwise |
| `nInterface` | the interface-cell count |

They match the definitions of leia's `leiaTestDavofNormal` (`E_L1_N`, `E_L2_N`,
`E_LINF_N`, `E_AREA_REL`), so `getResults.py`'s `results_sphereNormal3D.csv`
can be merged into the DAVOF order table with
`leia/workflow/scripts/make_davof_normal_table.py --tpf-results ...`.

Run (sourced OpenFOAM-v2606 and this clone's `etc/davof-env.sh`):

```
./runAll.sh                          # genCases.py -> 24 cases under Cases/ -> getResults.py
python3 genCases.py --N 20 40 --schemes plicRDF        # a subset
pytest run/benchmark/reconstruction/sphereNormal3D     # oftest ladder, plicRDF + gradAlpha, N <= 80
pytest run/benchmark/reconstruction/sphereNormal3D -m slow
```

`system/simulationParameter` carries the axes (`N`, `RECONSCHEME`,
`SETALPHAMETHOD` = `cutCellIso` | `cutCellImpFunc`); `genCases.py` and the pytest
rewrite it. `setAlphaMethod cutCellIso` initialises alpha from the linear
interpolant of the sphere function at the mesh points (the linearisation
comparable to leia's `linearInterpolant` DAVOF state); `cutCellImpFunc` cuts the
curved surface (bisection roots, refined cap), an alpha close to leia's exact
`exactSphere` state (closed-form face fractions, quadrature cell fractions).
Compare ORDERS across the two codes, values only within one `setAlphaMethod`;
the like-for-like value comparison is leia's in-app cross-check, which runs the
same plicRDF / gradAlpha / isoAlpha algorithms (OpenFOAM's libgeometricVoF) on
the identical alpha field.

Outputs: `results_sphereNormal3D.csv`, `orders_sphereNormal3D.csv` (least-squares
and pairwise orders per scheme x setAlphaMethod), `sphereNormal3D_L2.png`;
`results_pytest_sphereNormal3D.csv` from the pytest.
