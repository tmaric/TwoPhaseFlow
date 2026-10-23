# nonOrthogonalCavity2D

Free standing mode in a closed rigid square cavity on a smoothly distorted,
non-orthogonal mesh, solved with `acousticWaveFoam`. No PML, no interface,
no forcing and no added damping, so the case isolates the spatial operator
on non-orthogonal faces.

## Configuration

- Air in the square `0 < x, y < CAVITY_L` (0.343 m), rigid walls.
- Exact mode `p = PA cos(m pi x/L) cos(n pi y/L) cos(omega t)` with
  `(MODE_M, MODE_N) = (2, 1)`, initialized at the cell centres; `pDot`, `U`
  and the auxiliary fields start from zero.
- `blockMesh` builds `CELLS x CELLS` cells; `prepareCase --distort` then
  moves the points by
  `x' = x + DISTORTION L/(2 pi) sin(pi x/L) sin(2 pi y/L)` and the analogous
  `y'`. The walls stay straight. `DISTORTION = 0.15, 0.3, 0.6` give a maximum
  non-orthogonality of 17, 33 and 62 degrees (average 7, 13 and 25) at
  `CELLS = 40` (see `log.checkMesh`).
- `LAPLACIAN` selects the Laplacian and surface-normal gradient schemes.
  Only `orthogonal` runs at present: the solver scales the Laplacian matrix
  by `rho`, and OpenFOAM refuses to scale a matrix that carries the explicit
  non-orthogonal correction of `corrected` ("cannot scale a matrix containing
  a faceFluxCorrection"), for all three time updates.

## Run

```bash
./Allrun
./runStudy            # A: refinement at four distortions; B: long time
```

`Allrun` writes the cell centres and volumes (`writeCellCentres`,
`writeCellVolumes`), runs the solver with probes at every cell, and calls
`postprocess_cavity2d.py`: E2 of pressure and velocity against the exact
mode, the modal amplitude growth per period and frequency error, and the
drift of the physical energy. `summarize_nonorth_study.py` writes
`studies/studyA.csv` and `studies/studyB.csv`.

## Non-orthogonal correction

- `LAPLACIAN=orthogonal` drops the tangential gradient on non-orthogonal
  faces and does not converge on the distorted meshes.
- `LAPLACIAN=corrected` adds the explicit correction of OpenFOAM's scheme.
  The monolithic update converges it within each step with
  `NONORTH_CORRECTORS` block solves that reuse the factorization; the
  segregated update uses its corrector loop. Second order in space and
  time, but the corrected Laplacian is not symmetric: its complex
  eigenvalues make the lossless wave equation grow slowly on strongly
  distorted meshes (`amplification_cavity.py`).
- `LAPLACIAN=uncorrected NONORTH_MODE=symmetric` (monolithic; serial or
  parallel with processor patches, `NPROCS=2` or `4`) lets
  the solver apply the symmetrized correction `(C + C^T)/2`, which keeps the
  operator self-adjoint; `NONORTH_MODE=explicit` applies `C` itself and
  reproduces `corrected`.
- `WALL_LAYERS=1 TRANSITION_FRACTION=0.15` keeps the wall cells orthogonal;
  skewed wall cells with a zero-gradient wall reduce the order to one.

Studies: `runCorrectedStudy` (lagged correction), `runWallLayerStudy`,
`runDefectCorrectionStudy` (converged correction), `runSymmetricStudy`;
`summarize_corrected_study.py` and `summarize_wall_layer_study.py` write the
tables; `amplification_cavity.py OUTDIR` forms the (p, q) amplification
matrix and its spectrum.
