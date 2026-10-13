# pmlPulse2D

A Gaussian pressure pulse in a homogeneous gas leaves a square physical
domain through Cartesian PMLs on all four sides, solved with
`acousticWaveFoam`. The case exposes the PML to all angles of incidence,
its edges and corners included, and follows the solution long after the
pulse has left, to compare the monolithic and segregated updates, with
zero added damping.

## Configuration

- Physical domain `|x|, |y| < HALF_WIDTH` (0.05 m), surrounded by PMLs of
  thickness `PML_THICKNESS` (0.025 m), cubic profile (`PML_ORDER` 3) with
  `SIGMA_MAX` (5.05e5 1/s, continuous normal-incidence round trip 1e-8),
  sound-hard outer boundary. Uniform cubic cells of size `CELL_SIZE`.
- `DIM=3` gives a cube with PMLs on six sides; `DROP_RADIUS > 0` adds a
  centred water cylinder (2D) or sphere (3D).
- `ENLARGED=1` removes the PML and extends the domain to
  `ENLARGED_HALF_WIDTH` (0.115 m), the reference without truncation. Its
  cell centres coincide with those of the PML grid when
  `(ENLARGED_HALF_WIDTH - HALF_WIDTH - PML_THICKNESS)/CELL_SIZE` is an
  integer.
- Initial pressure `P0 exp(-|x - x0|^2/(2 s^2))` with `s = PULSE_WIDTH`
  (4 mm) at `SOURCE` (15, 10) mm, set by `setExprFields`; all other fields
  start from zero.
- `COURANT = c dt/h`; `TIME_METHOD` is `blockCoupled` (monolithic) or
  `thetaNewmark` (segregated with local implicit PML blend).

## Run

Load OpenFOAM-v2606 and `scripts/bashrc`, build `acousticWaveFoam`, then

```bash
./Allrun
```

Every entry of `caseParams.sh` can be overridden, e.g.
`CELL_SIZE=2e-3 END_TIME=0.05 SAMPLE_INTERVAL=1e-4 ./Allrun` for a long run.

## Monitoring

`cellProbes` records `p`, `pDot` and `U` at every physical cell;
`domainProbes` (`DOMAIN_PROBES=1`) records `p`, `v` (psi) and `u` (w) at
every cell, PML included. `fieldMinMax` is not used: it loads the VoF
library a second time, registers the `isoAlpha` reconstruction scheme twice
and aborts the solver at exit.

## Post-processing

`postprocess_pml.py` writes `postProcessing/pmlVerification/`:

- `E2` of pressure and velocity against the free-space solution (2D
  Hankel transform, 3D closed form, pressure only), time-integrated over the
  samples up to 0.4 ms; skipped with a drop;
- `energy.csv`: the physical energy `E_ac` (p, U) and the discrete wave
  energy `E_h` (pDot, p) of the physical domain;
- `maxima.csv`: the maxima of `|p|`, `|psi|`, `|w|` over the whole domain;
- growth rates over the last quarter of the run, each classified as
  decaying, bounded or growing (factor outside `[1/1.05, 1.05]`).

## Amplification matrices (study C)

Both time updates are linear one-step maps of the cell states. With `p`,
`p_0`, `pDot`, `v` and `u` read from the start time, `amplification.py`
forms the matrix `G` column by column: one solver step from `e_j + r`
minus one step from `r`, for a fixed random state `r` (unit states can
break down the iterative pressure solve). The state is (p, pDot, psi, w)
per cell and, for the segregated update, whose PML pressure equation uses
three time levels, also p^{n-1}. A step from a random state checks
linearity. `analyse_amplification.py` reports the spectral radius, the
eigenvalues outside the unit circle, the Jordan structure at g = 1 on the
active state (w in cells without PML stays zero and is removed), the norms
of G^n up to n = 2^20, the semidiscrete spectrum of the monolithic update,
and the explicit-damping radius |1 - sigma dt| of the outermost PML cell.
`SEGREGATED_IMPLICIT_DAMPING=1` selects the diagnostic segregated variant
that damps w implicitly.

```bash
./runSpectralStudy     # 34 matrices on the 16 x 16 mesh, studies/C/studyC.csv
```
