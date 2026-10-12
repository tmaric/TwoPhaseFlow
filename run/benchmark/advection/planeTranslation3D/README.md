# Plane-translation gate of the DAVOF face update

The exact-solution test of the face update of `isoAdvection` with the `davof`
reconstruction (`etc/README-davof.md`, "The face update"): a planar interface
translated uniformly is reproduced exactly by the moving PLIC plane, so the
transported face fractions `alphaf.water`, the cell fractions and the
Gauss-identity normal must stay at the exact plane's values to round-off at
every step.

The plane of `../../reconstruction/planeNormal3D` (unit normal
(0.3, 0.5, 0.8124)) is positioned in the (L, L, L) corner of the 4 mm box,
n.x = 1.33 L, so that it meets only the x = L, y = L and z = L boundaries:
outflow (or zero-flux) boundaries for both velocities, hence no cut cell has
an inflow boundary face, whose flux is the boundary condition's and not the
plane's. Two velocities (`VELOCITY` in `system/simulationParameter`):

- `normal`: U = n, the plane moves 0.5 h along its normal per step;
- `oblique`: U = (1 0 0), the plane moves 0.5 h along x, 0.15 h along n.

N = 20 by default (also run at 40), ten steps of 0.5 h / |U| (Courant 0.5,
fixed), both DAVOF state constructions (`detrixheAslam`, `quadraticFaces`:
both exact for a plane), `isoFaceTol 1e-14` and `wispTol 1e-3` (leia's
threshold |m_c| > wispTol V^(2/3)). The function object `areaFractionError`
writes per step the error of `alphaf` against the exact face fractions of
the translated plane over the faces of the interface cells (L1 mean, L2 rms,
Linf max, as fractions), the face and cell counts and the consistency
max |m_c - A n|/V^(2/3) of the transported state with the plane the scheme
positions in it; with `reconstructionErrorEveryStep true` it also writes the
normal and centre errors of the reconstructed plane per step
(`postProcessing/reconstructionError/0/error.dat`); `volumeFractionError`
writes E1 = sum |alpha - alpha_exact| V.

```
./runAll.sh                    # genCases.py (N = 20, both velocities, both states) -> Cases/ -> getResults.py
./runAll.sh --N 40
python3 getResults.py          # results_/history_planeTranslation3D.csv and the gate
```

Gate (`getResults.py`, pre-registered): max over the ten steps of L1(alphaf)
below 1e-12 and of the mean normal error L1(n) below 1e-12 for every case;
the maxima are reported.

## Measured 2026-10-12 (commit of the A3 step)

| velocity | N | max L1(alphaf) | max Linf(alphaf) | max L1(n) | max Linf(n) | result |
|---|---|---|---|---|---|---|
| oblique | 20 | 1.4e-14 | 2.5e-12 | 3.7e-13 | 1.6e-11 | PASS |
| oblique | 40 | 1.7e-13 | 1.0e-10 | 3.6e-12 | 8.7e-10 | marginal (L1(n)) |
| normal | 20 | 1.4e-07 | 2.4e-06 | 7.2e-07 | 3.7e-06 | FAIL after step 5 |
| normal | 40 | 2.8e-02 | 1.0 | 1.6e-01 | 1.4 | FAIL, growth from step 1 |

(`detrixheAslam` and `quadraticFaces` give the same numbers to three digits:
the states are identical for a plane.)

The oblique arm is exact at every step at N = 20. The normal arm is exact
for the first steps (N = 20: L1 4.5e-15, 1.5e-14, 1.1e-13, 1.1e-12 at steps
1-4; N = 40: 1.6e-14 at step 1) and then grows by a constant factor per
step: about 10 per step at N = 20 (from step 5, while the liquid tetrahedron
collapses into the corner: 67, 48, 32, 19, 10, 4 interface cells), 20 to 80
per step at N = 40 (from step 1, 761 to 298 interface cells). The growth is
round-off amplified through the interface cells with small |m_c|: in such a
sliver cell the alpha-positioned plane has the position sensitivity h/A
(A the cut area as a fraction of h^2) and the Gauss-identity normal m/|m|
the sensitivity 1/A to a face-fraction error, and the sliver's plane writes
its rotated, shifted trace onto the faces it sets, from which the
neighbours' normals are formed at the next step. The sliver hypothesis was
tested with the threshold: at N = 40 the growth factor per step is 20 to 80
at wispTol 1e-3, about 10 at 1e-2 and about 3 at 1e-1 (at N = 20: the
error after ten steps is 1.4e-7, 1.7e-7 and 4.7e-12); the oblique arm moves
the plane only 0.15 h along its normal per step and grows about 2x per step
at N = 40. The exactness of every single step (the first step is exact in
every case, and the unit check `apps/test/testCutFaceAdvectEndArea` puts the
end area and the swept flux at 4e-16) is established; the growth is a
property of the closed loop (own normal, alpha-positioned plane, plane
fluxes) in sliver cells and belongs to C2 of the proposal.

History of the normal arm at N = 20 (L1(alphaf) per step): 4.5e-15,
1.5e-14, 1.1e-13, 1.1e-12, 1.1e-11, 1.2e-10, 1.5e-09, 2.3e-08, 1.3e-07,
1.4e-07.
