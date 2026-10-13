# Sphere-translation advection benchmark

A sphere of radius 0.15 at (0.5, 0.5, 0.5) in the 2 x 1 x 1 box is carried
by the uniform velocity (1, 0, 0) to x = 1.5 in the end time 1 (one diameter
of travel every 0.3; inlet at x = 0 with alpha 0, outlet at x = 2, the lateral
boundaries carry no flux). N cells across the unit height (2N along x,
N = 32, 64, 128; the time step Co h / |U|, `genCases.py --Co`, default
0.2 since 2026-10-13, the cluster CSVs in this directory at 0.5, fixed;
Co is written into the case and the results CSV). Rows
(`genCases.py`): `isoAdvection` with the `davof` reconstruction on the DAVOF
state of `initDavofFields` (`davof-quadraticFaces`, `davof-detrixheAslam`),
and with `plicRDF` and `isoAlpha` on the `setAlphaField` state. Function
objects: `areaFractionError` every step (the error of the transported face
fractions `alphaf.water` against the exact face fractions of the translated
sphere, `davofMethod exactSphere`, over the faces of the interface cells of
either state; the consistency max |m_c - A n|/V^(2/3); the reconstruction
error of the scheme against the translated sphere at the end time, normal
and centre), `volumeFractionError` every step (E1 = sum |alpha -
alpha_exact| V against the translated sphere, the mass error, the
boundedness error). The runs end at the end time whatever the state does.

```
./runAll.sh --N 32 64          # the local ladder (N = 128: 4.2e6 cells)
./runAll.sh --jobs 24          # the cluster list, 24 cases at a time (../run-benchmark.sbatch)
./runAll.sh --case Cases/davof-quadraticFaces-exact/N32   # one case, no generation
python3 genCases.py --list     # the case list without generating
python3 genCases.py --Co 0.5   # the time step of the cluster runs of 2026-10-13
python3 getResults.py          # results_/orders_/history_/growth_sphereTranslation3D.csv
python3 getResults.py --cases DIR --suffix _laptopCo0.2   # another case tree into sibling CSV names
```

Rows since 2026-10-13 (`../advectionRows.py`; a row is scheme, davofMethod,
normalModel, applyTo, sliverTol, maxAngle in `system/simulationParameter`,
the davof scheme's normal model of `etc/README-davof.md`):
`davof-<method>` (every interface cell the DAVOF normal, `cellLocal`),
`davof-<method>-Youngs<c>` (the mixed normal after Basilisk's
Youngs-centred scheme: slivers |m_c| < c h^2 take `fvc::grad(alpha)`;
`youngs`, check, maxAngle 180), `davof-<method>-exact` (control: every
interface cell the exact normal of the translated sphere at the plane
centroid, alphaf still transported by the moved plane),
`davof-<method>-plicRDF` (control: every normal from plicRDF), optionally
`davof-<method>-sliverPlicRDF<c>`, the ring rows `ringArea` (the
area-weighted ring normal in every interface cell), `ringArea-check` (in
the slivers |m_c| < 0.1 h^2 and the cells deviating by more than 15
degrees), `ringFit` and `ringFit-check` (the ring quadric fit), and the
comparators `plicRDF` and `isoAlpha`. The cluster list is the state
quadraticFaces, c = 0.01, 0.1, 1, both controls, the four ring rows and
the comparators at N = 32, 64, 128 (36 cases, one node, `sbatch
--export=ALL,REPO=...,BENCH=$PWD ../run-benchmark.sbatch`);
`--davofMethod detrixheAslam` gives the same rows on the Detrixhe-Aslam
state.

## Measured 2026-10-13 on the laptop at N = 32, Co 0.2 (`*_laptopCo0.2.csv`)

After the normal models, the phantom-plane fix of `davof::positionPlane`
(the cause of the unbounded cluster runs of the same day: a face-fraction
interface cell with alpha at 0 or 1 inherited the previously positioned
cell's plane) and the two flux bounds (`etc/README-davof.md`, which has
the full table): every row is bounded (min alpha -9.8e-9 for the own
normal, round-off otherwise). End values at t = 1 (L1(alphaf), the normal
error L1(n), the centre error L1(x), E1, interface cells):
`davof-quadraticFaces` 2.08e-1, 1.29, 1.54e-1, 5.93e-3, 4910 (the own
normal still loses the state, L1(alphaf) 4.4e-3, 5.5e-3, 8.1e-3, 1.6e-2,
6.1e-2 at steps 1, 2, 3, 5, 10); `davof-quadraticFaces-exact` 1.04e-2,
3.64e-3, 9.83e-4, 2.76e-4, 426 (4.4e-3, 7.6e-3, 5.7e-3, 6.0e-3, 6.7e-3);
`ringArea` 1.78e-2, 4.06e-2, 1.35e-3, 3.85e-4, 424 (5.2e-3, 6.1e-3,
6.0e-3, 7.0e-3, 7.7e-3: the loop closed, the interface-cell count
constant); `ringArea-check` 4.49e-2, 1.32e-1, 3.36e-3, 8.13e-4, 487;
`ringFit` 7.93e-2, 2.05e-1, 1.32e-2, 1.15e-3, 687 (matches ringArea to
step 80, then 176 fit fallbacks); `ringFit-check` 9.90e-2, 3.01e-1,
1.69e-2, 1.23e-3, 765; `plicRDF` -, 3.92e-2, 2.05e-3, 3.18e-4; `isoAlpha`
-, 5.94e-2, 2.17e-3, 3.55e-4. The ladder (N = 64, 128) and the orders
are for the cluster. `growth_sphereTranslation3D.csv` holds
per davof case L1(alphaf) of the first 20 steps and the growth factors
L1(k)/L1(k-1); `history_...csv` the sliver count of the advection step's
reconstruction per step. If the `exact` row is stable and first order the
face update is sound and the own-normal feedback is the culprit (the plane
gate says so at N = 20, `../planeTranslation3D/README.md`).

`getResults.py` prints per row and N the final L1(alphaf), the normal error
(L1, L2), the centre error (L1, L2), E1 with the pairwise orders over N, the
mass error and the alpha bounds of the last step; `history_...csv` carries
the per-step `areaFractionError` rows of the davof cases.

Prediction of the proposal (finding (v), the iterated regeneration): the
quadratic initial state loses its second order within the first steps (the
plane clipping is q = 1), so L1(alphaf) after one diameter of travel is
first order, the normal first order, the centre second order.

## Measured 2026-10-12 (commit of the A4(b) step), N = 32 and 64

End of the run (t = 1, 3.3 diameters of travel); the reconstruction errors
are the mean |n - n_exact| and the mean distance of the reconstructed centre
from the translated sphere over the interface cells, E1 the volume error,
the orders pairwise 32 -> 64:

| row | N | L1(alphaf) | L1(n) | L1(x) | E1 | Emass | min alpha |
|---|---|---|---|---|---|---|---|
| plicRDF | 32 | - | 2.50e-2 | 1.67e-3 | 4.93e-5 | 5.6e-17 | -4e-15 |
| plicRDF | 64 | - | 1.43e-2 (0.81) | 4.11e-4 (2.02) | 1.76e-5 (1.49) | 2.9e-17 | -2e-15 |
| isoAlpha | 32 | - | 4.36e-2 | 1.67e-3 | 4.45e-5 | 2.3e-17 | -9e-15 |
| isoAlpha | 64 | - | 2.80e-2 (0.64) | 4.45e-4 (1.91) | 1.34e-5 (1.73) | 2.4e-17 | -2e-15 |
| davof-quadraticFaces | 32 | 1.42e-1 | 1.39 | 5.36e-1 | 8.44e-2 | 2.7e-3 | -2.0 |
| davof-quadraticFaces | 64 | 1.25e-1 (0.18) | 1.38 (0.01) | 5.69e-1 | 8.16e-2 (0.05) | 4.8e-3 | -2.2 |
| davof-detrixheAslam | 32 | 1.36e-1 | 1.40 | 5.33e-1 | 8.51e-2 | 2.4e-3 | -1.5 |
| davof-detrixheAslam | 64 | 1.24e-1 (0.14) | 1.38 (0.01) | 5.68e-1 | 8.12e-2 (0.07) | 4.6e-3 | -2.5 |

The comparators behave as expected (normal first order, centre second
order, E1 between first and second order, mass conserved, bounded). The
davof loop is unstable: at N = 32 (R/h = 4.8) the face-fraction error of
the quadratic state is 1.2e-2 after one step (Linf 0.12), 2.6e-2, 3.5e-2
(Linf 1), 8.4e-2, 0.20, 0.27, 0.32, 0.37, 0.39, 0.40 at steps 1 to 10,
alpha falls below -1e-6 at step 4 (min -0.5, `clip false`), the
interface-cell count grows from 452 to 9461 at step 10 and 58 000 of
65 000 cells at the end, and the consistency max |m_c - A n|/V^(2/3) is
about 1 from the first step and 3 afterwards; the Detrixhe-Aslam state
does the same (1.1e-2, 3.7e-2, 5.5e-2, 9.6e-2, 0.20, ... at steps 1 to
10). At N = 64 the picture is the same (final L1 0.125, "order" 0.18). The
per-step histories are in `history_sphereTranslation3D.csv`.

The mechanism is the one measured in the plane gate
(`../planeTranslation3D`): interface cells with small |m_c| (slivers just
above wispTol) have an ill-conditioned alpha-positioned plane and an
ill-conditioned own normal m/|m|, and their rotated planes write onto the
faces they set; on the sphere the seed is the O(h/R) plane-clipping error
of every step instead of round-off, and the amplification (10 to 80 per
step in the gate at wispTol 1e-3) leaves no first-order plateau to measure.
The rules of the face update were varied during the day without changing
the first three steps (the second ring setting the fluxes always, only
where the plane crosses the face, or never): the breakdown comes from the
surface cells' own planes. The prediction of the proposal (first-order
L1(alphaf) after one diameter, first-order normal, second-order centre)
is therefore neither confirmed nor refuted by this run; the loop needs
the stabilisation the proposal lists (Stage 1: the clipping surface that
reproduces the cell's area vector, the projection) or a threshold for
plane-setting cells above the interface-cell threshold before the orders
can be measured.

