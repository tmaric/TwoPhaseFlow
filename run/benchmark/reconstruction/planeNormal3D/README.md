# Plane reconstruction gate

The exact-solution gate of the `davof` reconstruction scheme: a tilted plane
(normal (0.3, 0.5, 0.8124), position (0.00213, 0.00171, 0.00197) m, the same
plane as leia's `cases/davof/planeNormal3D`) through a 4 mm box of N^3 cells,
N = 20 by default. Every linear interpolant is exact for a plane, so both DAVOF
states (`davof-detrixheAslam`, `davof-quadraticFaces`) carry the exact face and
cell fractions, and the Gauss-identity normal, the alpha-positioned plane and
the explicit plane position must come out at round-off:

- `LNormalDiffInf` (max |n - n_exact| over the interface cells) below 1e-12,
- `LCentreDiff2`/h (rms distance of the reconstructed centre from the plane,
  relative to h) below 1e-12 and `LCentreDiffInf`/h below 1e-10 (the case's
  `isoFaceTol` is 1e-14 for that; the maximum is the floor of the
  alpha-positioned plane in corner-sliver cells: measured 2026-10-12 at
  N = 20, three cells with 1 - alpha = 1.8e-6 and a cut polygon of 5e-4 h^2
  set 1.3e-11 h, the rms is 7.7e-13 h and the median 1e-13 h; leia's explicit
  position has the same 1e-11 h floor),
- the explicit plane's defects printed by the scheme, max |alpha_plane - alpha|
  and max |A n - m|/V^(2/3), below 1e-12.

`getResults.py` checks the five numbers for every `davof` row and exits
non-zero on a failure. Measured 2026-10-12 (N = 20, both states):
|n - n_exact| max 8.4e-13 / 7.7e-13, rms distance 7.7e-13 h, max 1.3e-11 h,
plane defects 4.6e-15 / 2.7e-13; plicRDF on the same alpha: normal rms
1.4e-4, centre rms 2.2e-12 m. `plicRDF`, `gradAlpha` and `isoAlpha` run on the same
alpha (the `davof-detrixheAslam` state and the case's own `cutCellIso`) and are
reported, not gated.

Same scripts as `../sphereNormal3D` (`genCases.py`, `runAll.sh`,
`getResults.py`, keyed by the directory name); `system/simulationParameter`
carries the axes (`N`, `RECONSCHEME`, `SETALPHAMETHOD`, `DAVOFMETHOD`).

```
./runAll.sh                      # N = 20, davof on both states, comparators
python3 getResults.py            # results_planeNormal3D.csv, orders_planeNormal3D.csv, the gate
```
