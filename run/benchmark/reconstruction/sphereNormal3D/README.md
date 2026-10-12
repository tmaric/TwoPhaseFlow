# Sphere-normal reconstruction ladder

Static convergence of the PLIC interface normal reconstructed by `davof`,
`plicRDF`, `gradAlpha` and `isoAlpha` from the volume fraction of a sphere
(R = 1 mm in a 4R box, centre off every mesh symmetry plane) on uniform hex
meshes with N = 20, 40, 80, 160 cells per side (R/h = 5, 10, 20, 40). One
reconstruction step per case (`reconstructInterface`, `stopAt writeNow`).

The `davof` scheme (`src/VoF/reconstructionSchemes/plicSchemes/davof`)
recovers the normal from the face fractions `alphaf` of the cell's own faces
(the Gauss identity m_c = sum_f alphaf S_f^out) and positions the plane by
alpha; it needs the DAVOF state, which `initDavofFields` writes
(`SETALPHAMETHOD davof`, `DAVOFMETHOD detrixheAslam | quadraticFaces |
exactSphere` in `system/simulationParameter`; the Allrun runs it before
`reconstructInterface`, whose `setAlphaMethod davof` reads alpha and alphaf
from `0/`). The comparators run on the same alpha (genCases.py skips them on
`davof-quadraticFaces`, whose cell fractions equal `davof-detrixheAslam`'s).
The explicit plane position and its realizability defects are written by the
scheme (`pPlane`, `planeVolDefect`, `planeAreaDefect`) and their maxima
collected from the log (`maxPlaneVolDefect`, `maxPlaneAreaDefect`).

The error columns were appended to `reconstructionError` (`src/postProcessing`)
for this benchmark, after the existing twelve, which are unchanged:

| column | meaning |
|---|---|
| `LNormalDiff1`, `LNormalDiff2`, `LNormalDiffInf` | mean, rms, max over the interface cells (`|normal| != 0`) of `|n - n_exact|`, `n` the unit PLIC normal, `n_exact` the sphere normal at the reconstructed centre (both point into alpha = 1) |
| `areaSum` | `sum |normal|` over the cells: `normal_` is the PLIC area vector, so this is the reconstructed interface area |
| `areaExact` | `4 pi R^2` for `type sphere`, 0 otherwise |
| `nInterface` | the interface-cell count |
| `LCentreDiff1`, `LCentreDiff2`, `LCentreDiffInf` | mean, rms, max over the interface cells of the distance of the reconstructed centre from the surface (exact `|R - |x - o||` for the sphere, the Newton-projected distance to the implicit function otherwise; since 2026-10-12) |

They match the definitions of leia's `leiaTestDavofNormal` (`E_L1_N`, `E_L2_N`,
`E_LINF_N`, `E_AREA_REL`), so `getResults.py`'s `results_sphereNormal3D.csv`
can be merged into the DAVOF order table with
`leia/workflow/scripts/make_davof_normal_table.py --tpf-results ...`.

Run (sourced OpenFOAM-v2606 and this clone's `etc/davof-env.sh`):

```
./runAll.sh                          # genCases.py -> the default ladder under Cases/ -> getResults.py
./runAll.sh --N 20 40 80 --setAlpha davof-detrixheAslam davof-quadraticFaces davof-exactSphere
python3 genCases.py --N 20 40 --schemes plicRDF --setAlpha cutCellIso   # a subset
pytest run/benchmark/reconstruction/sphereNormal3D     # oftest ladder, plicRDF + gradAlpha, N <= 80
pytest run/benchmark/reconstruction/sphereNormal3D -m slow
```

`system/simulationParameter` carries the axes (`N`, `RECONSCHEME`,
`SETALPHAMETHOD` = `cutCellIso` | `cutCellImpFunc` | `davof`, `DAVOFMETHOD`);
`genCases.py` and the pytest rewrite it. `setAlphaMethod cutCellIso`
initialises alpha from the linear interpolant of the sphere function at the
mesh points (the linearisation comparable to the `davof-detrixheAslam` state);
`cutCellImpFunc` cuts the curved surface (bisection roots, refined cap), an
alpha close to the exact `davof-exactSphere` state (closed-form face fractions,
quadrature cell fractions). Values are comparable only within one
`setAlphaMethod`; across the two codes compare the orders, and the values of
the `davof-*` states, which are the same constructions as leia's
`detrixheAslam`, `quadraticFaces` and `exactSphere` sources.

Outputs: `results_sphereNormal3D.csv`, `orders_sphereNormal3D.csv` (least-squares
and pairwise orders per scheme x setAlphaMethod, normal and position),
`sphereNormal3D_L2.png`; `results_pytest_sphereNormal3D.csv` from the pytest.

Measured 2026-10-12 (N = 20, 40, 80; least-squares orders of the normal error,
L1 / L2, and of the position error L2, with the L2 normal error at N = 80):
see the commit message of the A4(a) step and `orders_sphereNormal3D.csv`.
Those CSVs were measured with `WISPTOL 1e-6`; since the A3 step the
parameter files (here and in `../planeNormal3D`) carry leia's `1e-3`, which
only drops the sliver cells from the interface-cell count (the orders and
the L2 values of the sphere ladder change in the last digit at most; the
ladder has not been re-run).
