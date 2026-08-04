# homogeneousPlaneWave2D

Source-free analytical Helmholtz verification in a homogeneous two-dimensional
domain. The case supports oblique exact Dirichlet data and an aligned wave with
exact mixed Dirichlet--Neumann data on orthogonal and smoothly warped meshes.

Run one case with `./Allrun`. Override controls through the environment, for
example:

```bash
CELLS_PER_WAVELENGTH=48 MESH_FAMILY=warped BOUNDARY_MODE=mixed ./Allrun
```

Run the complete convergence study with `./run_convergence.sh`. Results are
written to `postProcessing/convergence`.

