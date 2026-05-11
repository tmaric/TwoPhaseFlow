# WedgeLevitatorGorkovHole

This case combines:

- the `WedgeLevitatorDropHole` cfMesh hard-boundary drop mesh, where the drop
  is represented by a `dropWall` patch;
- the `WedgeLevitatorGorkovAlpha` radiation-force postprocessing layout.

The `dropRadiationForce` function object integrates `pr` directly on
`dropWall` and writes `postProcessing/dropRadiationForce/force.dat` with the
same column layout used by `GorkovAlpha`.

The radius sweep writes `solidBallGorkovStudy/summary.tsv` and
`solidBallGorkovStudy/acousticForceVsRadius.png`.

Run:

```sh
./Allrun
```

Optional radius sweep:

```sh
./runSolidBallGorkovStudy
```
