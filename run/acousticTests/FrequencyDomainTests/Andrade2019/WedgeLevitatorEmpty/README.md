# WedgeLevitatorEmpty

This case combines:

- the same `gmsh`/`gmshToFoam` outer levitator mesh workflow used by
  `DropAlpha`;
- no drop/hole inside the acoustic domain;
- a coded function object that writes the Gorkov potential and acoustic
  radiation force fields.

The function object writes `EpGorkov`, `EkGorkov`, `GorkovPotential`, and
`acousticRadiationForce` at the latest time. The force is
`-grad(GorkovPotential)` and the potential uses the probe particle properties
in `caseParams.sh`.

It also writes `postProcessing/gorkovPointForce/force.dat`, where the force is
sampled from the nearest cell to `(AXIS_RADIUS, PARTICLE_CENTER_Y, 0)`.

Run:

```sh
./Allrun
```

To compare the empty-field Gorkov point-force estimate with the hard-hole
integrated force from `WedgeLevitatorGorkovHole`, run:

```sh
./runEmptyGorkovStudy
```

This writes `emptyGorkovStudy/summary.tsv`,
`emptyGorkovStudy/gorkovVsHoleForce.tsv`, and
`emptyGorkovStudy/gorkovVsHoleForce.png`.
