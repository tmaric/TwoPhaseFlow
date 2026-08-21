# homogeneousPlaneWave2D

Source-free analytical Helmholtz verification in a homogeneous two-dimensional
domain. The case supports oblique exact Dirichlet data and an aligned wave with
exact mixed Dirichlet--Neumann data on orthogonal meshes and meshes with a
non-orthogonal interior.  The `warpedInterior` family keeps the requested
number of boundary-cell layers orthogonal and smoothly ramps the deformation
towards the interior.

Run one case with `./Allrun`. Override controls through the environment, for
example:

```bash
CELLS_PER_WAVELENGTH=48 MESH_FAMILY=warpedInterior \
    ORTHOGONAL_BOUNDARY_LAYERS=1 BOUNDARY_TRANSITION_FRACTION=0.15 \
    BOUNDARY_MODE=mixed ./Allrun
```

Run the complete convergence study with `./run_convergence.sh`. Results are
written to `studyResults/convergence`.
