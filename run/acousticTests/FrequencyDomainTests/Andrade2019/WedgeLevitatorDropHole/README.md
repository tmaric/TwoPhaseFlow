# WedgeLevitatorDropHole

This case combines:

- the `WedgeLevitatorEquilibriumHole` cfMesh hard-boundary drop mesh, where the drop
  is represented by a `dropWall` patch;
- the `WedgeLevitatorGorkovAlpha` radiation-force postprocessing layout.

The `dropRadiationForce` function object integrates `pr` directly on
`dropWall` and writes `postProcessing/dropRadiationForce/force.dat` with the
same column layout used by `GorkovAlpha`.

The hard-ball radius sweep writes `solidBallGorkovStudy/summary.tsv` and
`solidBallGorkovStudy/acousticForceVsRadius.png`.

The Gorkov comparison is now self-contained in this case folder. It first uses
the `emptyGmshTemplate` gmsh workflow to solve the empty levitator field and
postprocess the Gorkov point force for each radius. It then uses this DropHole
cfMesh workflow for the hard-ball acoustic force and compares both summaries
with `plotGorkovForceComparison.py`.

Run:

```sh
./Allrun
```

Optional radius sweep:

```sh
./runSolidBallGorkovStudy
```

Full empty-field Gorkov vs hard-ball comparison:

```sh
./runGorkovComparison
```

or for an explicit radius range:

```sh
./runGorkovComparison 20 80 5
```
