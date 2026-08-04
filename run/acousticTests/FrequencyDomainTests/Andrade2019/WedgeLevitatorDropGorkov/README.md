# Rigid-sphere Gorkov verification

This case compares the radiation force obtained by integrating the solver's
time-averaged radiation-stress traction over a resolved sound-hard sphere with
the Gorkov force for a rigid Rayleigh particle in a one-dimensional standing
wave. The traction contains the radiation-pressure and momentum-flux terms;
the latter has zero normal contribution on the sound-hard sphere to numerical
precision.

The axisymmetric cavity has equal transducer and reflector radii. The
reflector, cylindrical side wall, and sphere are sound hard. At the
transducer, the analytical pressure amplitude

`P(y) = Pa cos(k y)`

is imposed as a fixed complex-pressure boundary value. `Pa` therefore sets
the excitation; the field around the sphere is solved and is not prescribed.
The PML damping is zero for this closed-cavity verification.

Run one case with:

```sh
./Allrun
```

Run the radius study and generate its table and figure with:

```sh
./runGorkovStudy
```

Run the three-level near-sphere mesh study and the particle-position sweep with:

```sh
./runMeshConvergence
./runPositionSweep
```

The main outputs are `gorkovStudy/gorkov_comparison.tsv` and
`gorkovStudy/gorkov_comparison.png`. The refinement and position studies write
the corresponding consolidated files under `meshConvergence/` and
`positionSweep/`.

The default sweep covers radii from 30 to 150 micrometres
(`0.0139 <= ka <= 0.0694`). The near-sphere target cell size is one sixteenth
of the radius, with 128 segments on the generating semicircle.
Set `KEEP_CASES=1` when running the study to retain each generated mesh and
solution; otherwise the script keeps only the consolidated outputs.

Because the reference excitation of 1 Pa produces femtonewton forces, the lobe
and first zero-crossing cases were repeated with `PRESSURE_AMPLITUDE=1000`.
After division by the required factor of $1000^2$, the forces agree with the
1 Pa results to relative differences of $1.72e-8$ and $4.74e-9$. The retained
values are in `positionSweep/amplitude_scaling.tsv`.
