# The `davof` branch

Reference build of TwoPhaseFlow for the DAVOF (dual area/volume-of-fluid) benchmarks, which
live in the sibling `leia` clone (`<workspace>/leia`, branch `development`).

## What it is

- Base: the fork branch `feature/parasitic-currents` (5aa7e9a), i.e. upstream `de9826f` plus
  the `Foam::twoPhaseFlow` inline-namespace fix for the run-time-selection collision between
  this repository's `libVoF` and OpenFOAM's built-in `libgeometricVoF` (same class names in
  the same namespaces; libgeometricVoF is co-loaded whenever a case uses
  `libs (fieldFunctionObjects)`), and the `run/parasiticCurrents/` reproducer cases. Root
  cause and evidence: `run/parasiticCurrents/README.md`.
- Merged with upstream `master` `ba10702` (2026-09-08): the OpenFOAM-v2606 port and the
  heightFunction parallel fix.
- Case-syntax change inherited from the fix: the sampled surface is selected with
  `type reconstructedInterface;` (OpenFOAM's geometricVoF owns the key `interface`).

## Build

```
source $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc
export WM_NCOMPPROCS=24
./Allwmake
```

`Allwmake` sources `etc/davof-env.sh`, which redirects `FOAM_USER_APPBIN` / `FOAM_USER_LIBBIN`
(the only install variables the `Make/files` use) to `platforms/<WM_OPTIONS>/{bin,lib}` inside
this clone, and removes the account-wide `$HOME/OpenFOAM/$USER-<version>/platforms` from
`PATH` and `LD_LIBRARY_PATH`. No `Make/files` or `Make/options` is changed, so merging
upstream stays trivial. Nothing is written to, or deleted from, the account default.

`platforms/<WM_OPTIONS>/.openfoam-version` records the OpenFOAM version of the build.
v2512 and v2606 share `WM_OPTIONS`, so `davof-env.sh` refuses to run when the sourced
OpenFOAM differs from the stamp: `./Allwclean && ./Allwmake` under the version you want.

To run a binary of this clone in a fresh shell:

```
source $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc && . ./etc/davof-env.sh
```

## Verification of the collision fix

```
. etc/davof-env.sh; L=$FOAM_USER_LIBBIN/libVoF.so; G=$FOAM_LIBBIN/libgeometricVoF.so
nm -DC --defined-only $L | grep -c 'Foam::twoPhaseFlow::reconstructionSchemes'   # > 0
nm -DC --defined-only $L | grep -c ' Foam::reconstructionSchemes'                # 0
strong() { nm -DC --defined-only "$1" | awk '$2 ~ /^[TDBR]$/ { $1=""; $2=""; sub(/^  /,""); print }' | sort -u; }
comm -12 <(strong $L) <(strong $G) | wc -l                                        # 0
(cd run/parasiticCurrents/stationaryDroplet2D && ./Allclean && ./Allrun)          # "interFlow exit code: 0"
```

Known, accepted leftovers: `libVoFphaseModelCompressibleTurbulenceModels` still shares
template instantiations with OpenFOAM's `libVoFphaseCompressibleTurbulenceModels` (only
when loaded explicitly); `dynamicAlphaContactAngleJiang` is defined in two of this
repository's own libraries.

## The face update (isoAdvection with the davof reconstruction)

`isoAdvection` stores, when the reconstruction scheme is `davof` (or
`storeAreaFractions true;` in the alpha solver dictionary), the submerged
area of every face at the END of the step as the face fraction
`alphaf.<group>` (the registered field the `davof` scheme reads; created
upwind when none is registered). `cutFaceAdvect::timeIntegratedFaceFlux`
keeps the liquid area of the face at t = dt of its sweep
(`subFaceAreaEnd()`: the face fully on the side the plane has passed or
not reached, or the submerged area at the last vertex time plus the last
swept quadrilateral `initialArea + sign(Un0) (alpha + beta)` at dt; summed
over the fan of triangles of a non-planar face and scaled to the face area;
the stationary-plane branch stores `|subFaceArea|`; the `nShifts` warning
branch stores -1) and whether the plane cut the face at some time in
[0, dt] (`faceSwept()`). Unit check: `apps/test/testCutFaceAdvectEndArea`
(a unit cube's six faces, six sweep directions, both signs of Un0, five
offsets, nine dt: the end area and the swept flux against the exact clipped
area and its piecewise-quadratic time integral, both at 4e-16).

Rules of `isoAdvection::timeIntegratedFlux` (counts in one Info line per
step):

1. Before the sweep every face holds the upwind cell's alpha, boundary faces
   the owner's (`defaultAreaFractions`); after the bounding the faces no
   plane has set are re-evaluated with the new alpha.
2. The downwind faces of a surface cell: the cell's plane sets `alphaf` and,
   as before, `dVf`.
3. An upwind face of a surface cell: the upwind cell's plane sets it when
   that cell is a surface cell (its downwind face, rule 2); otherwise this
   cell's plane sweeps it and sets `alphaf`, and `dVf` where the plane can
   speak for the upwind cell W (`planeSpeaks`): W a partial cell below the
   interface-cell threshold (a sliver, whose upwind value phi alpha_W dt is
   neither the side value nor the swept part), or W pure and the plane
   crossing the face during the step (the swept part); a pure W the plane
   does not reach keeps the upwind flux (the side value of an extended
   plane can contradict a pure cell's content on a curved interface).
4. Boundary faces of a surface cell: the owner's plane sets `alphaf`; the
   outflow flux is the plane's as before, an inflow face keeps the boundary
   condition's flux.
5. The second ring: every face no plane has set whose cell shares a face
   with a surface cell is swept by the neighbouring plane whose
   end-of-step position is closest to the face centre (one and the same
   plane for a planar interface), which sets `alphaf` (the swept part or
   the side value) and `dVf` under the same `planeSpeaks` condition as
   rule 3 (inflow boundary faces excepted). A face the moving plane
   crosses during the step without belonging to a surface cell exists
   whenever a face-normal component of n is below the Courant number (the
   plane passes within Co h of a vertex ahead of both cells); without
   this rule the plane gate fails at the second step (the flux through
   such a face is 0 upwind instead of the swept part), and without the
   side value on the faces of a cell the plane enters it fails as well
   (the upwind default re-evaluated with the new alpha hands the cell's
   fraction to all its faces).
6. Processor patches: `alphaf` is exchanged with `dVf` (the upwind side's
   value, no sign change); the second ring does not cross processor
   boundaries.
7. A sweep that cannot set the face (`subFaceAreaEnd() < 0`) keeps the
   default and is counted as a fallback; values outside [0, 1] by more
   than 1e-12 are counted as clips and clamped.
8. The content bound (2026-10-13, `boundCellContent`, key `contentBound`,
   default true): for every non-surface cell W with faces whose flux a
   plane set (rules 3 and 5, `planeFlux`), the liquid leaving W through
   its faces may not exceed alpha_W V_W plus the liquid entering it in
   the step, nor the gas (1 - alpha_W) V_W plus the entering gas; when
   violated, the plane-set fluxes of W (its outflow faces) are moved
   towards their upwind values by one factor lambda in [0, 1] chosen so
   that the bound holds with equality (the upwind fluxes satisfy it for
   cell Courant numbers below one); the bounded faces and cells are
   counted in the Info line ("content bound on F faces of C cells").
   Tolerance 1e-12 V_W, so that an exactly translating plane never
   triggers it (the plane rows at N = 20 and 40, Co 0.2 and 0.5: 0 faces).
9. The conservative bounding of `limitFluxes` is extended by a second
   layer of cells when the face fractions are stored (key
   `boundSecondLayer`, default true): the second-ring sweep sets plane
   fluxes on faces whose other cell lies beyond the first layer
   (surface cells plus one), which `boundFlux` otherwise never visits;
   the cells beyond the first layer it corrects are counted ("bounding
   extended by a second layer: K cells corrected beyond the first
   layer"). `debugCell <label>` with `DebugSwitches { isoAdvection 1; }`
   prints the content bound of one cell.

The unbounded sphere of the cluster runs of 2026-10-13 (every davof row,
the exact-normal control included, min(alpha) at -Co within a few steps)
was NOT a flux rule but a reconstruction defect: `davof::positionPlane`
read `surfaceIteratorPLIC::cellStatus()` after `vofCutCell`, which
returns -1 or 1 for alpha below `isoFaceTol` or above 1 - `isoFaceTol`
BEFORE it cuts anything, so `cellStatus()`, `surfaceArea()` and
`surfaceCentre()` still held the previously positioned cell's plane. A
davof interface cell is selected by its face fractions and can have
alpha at 0 or 1 (a cell the interface has just left or entered); such a
cell became a phantom surface cell carrying its predecessor's plane,
which isoAdvection then fluxed (a full face of liquid out of an empty
cell 1.5 h outside the sphere, found with `debugCell` on cell 42190 of
the N = 32 case at step 6). Since 2026-10-13 `positionPlane` uses the
returned status (plicRDF never saw this because it positions planes
only in cells with surfCellTol < alpha < 1 - surfCellTol). Measured on
the exact-normal row of the sphere at N = 32, Co 0.2 (160 steps): before
the fix min(alpha) reaches -0.196 at step 5 and -1.83 by t = 0.5, with
L1(alphaf) 2.8e-2 at step 4 and 0.13 at step 5; the second layer alone
delays the event by one step (min -2e-25 at step 5, -0.095 at step 6,
-0.2 at step 7, 4714 corrections per step at the end); the content bound
on top keeps step 5 at round-off and still loses step 6 (-0.095, a cell
inside the first layer that `boundFlux` cannot repair in three passes
because its inflow came from a phantom plane). After the fix the row is
bounded with both bounds OFF (min -7.8e-15 over 160 steps, L1(alphaf)
4.4e-3, 7.6e-3, 5.7e-3, 6.0e-3, 6.0e-3 at steps 1-5, 1.04e-2 at t = 1,
the normal error L1 1.8e-5 at the end) and with both on (min -8.7e-15,
the content bound active on 19, 8, 25 faces at steps 2-4 and on none
after step 9; L1(alphaf) identical to four digits). The plane rows at
N = 20, Co 0.5, are bit-identical with and without the fix.

Nothing changes for the other reconstruction schemes (plicRDF after ten
steps of the deformation case is bit-identical with the build before this
change, alpha and interfaceNormal). Dictionary keys of the alpha solver:
`storeAreaFractions` (default: the scheme is davof), the davof scheme's
`wispTol` (default 1e-3 since this step, leia's value), `isoFaceTol`.

### Normal models (2026-10-13)

The normal along which the `davof` scheme positions the plane of an
interface cell is a run-time-selectable model (`src/VoF/davof/normalModels`,
base `davofNormalModel`, keyword `normalModel` of the alpha solver
dictionary, default `cellLocal`). `reconstruct()` forms m_c for every
interface cell, calls the model, which returns the normal it proposes per
interface cell, and positions the plane of every interface cell by alpha
along the selected normal (the model normal where it replaces m_c, m_c
otherwise). The models that read the planes positioned along m_c
(`needsPositionedPlanes()`: `exact`, which takes the exact normal at the
centroid of that plane, and `ringFit`, which fits the ring's planes) get
them in a first pass, and the cells whose normal changed are re-positioned
in a second one, as every model was until 2026-10-15; the other models
position once. The plane is the same either way (`vofCutCell` starts from
the cell's vertex values, nothing is carried over from an earlier cut), so
the results are bit-identical; the single pass saves one positioning per
interface cell for `cellLocal`, `youngs`, `plicRDF` and `ringArea`. Every
model takes `applyTo`:

- `all`: every interface cell with a model normal is replaced;
- `check`: only the slivers, |m_c| < `sliverTol` V_c^(2/3) (default 0.1),
  and the cells whose raw normal deviates from the model normal by more
  than `maxAngle` degrees (default 15, plicRDF's residual rule).

A cell wanted by the rule whose model normal is zero gets what the keyword
`noNormalFallback` says (base class, every model; since 2026-10-14):
`keep` (the default) keeps m_c, counted "kept m"; `ringArea` gives it the
ring normal of the cell, the sum of m_j over the point-neighbour ring
including the cell (the ringArea model's computation,
`davofNormalModel::ringAreaNormals()`, one implementation shared by the
model and the fallback; the exchange is collective and runs once per
reconstruction when the fallback is on), counted "ring fallbacks"; a cell
whose ring sum is zero as well keeps m_c. It matters for `plicRDF` (the
helper marks only the cells with alpha strictly inside (surfCellTol,
1 - surfCellTol), the davof scheme its interface cells by the face
fractions) and `youngs` (grad(alpha) can vanish); `exact` always has a
normal and the ring models' own normal is the ring sum. One Info line per
reconstruction: `davof: N interface cells (wispTol c), normalModel X
(applyTo a, sliverTol s[, maxAngle m][, noNormalFallback f]): S slivers, R
replaced[, K kept m][, F fit fallbacks][, G ring fallbacks], explicit plane
defects ...` (the keyword is printed only when it is not `keep`, so the
logs of every existing row are unchanged). `sliverFlag.<group>` is written
(1 in slivers). Models:

- `cellLocal`: the raw m_c/|m_c|, nothing replaced (default applyTo all).
- `youngs`: `fvc::grad(alpha)` (into alpha = 1 like m_c; default check).
- `plicRDF`: the normal of a helper plicRDF reconstruction of the same
  alpha, owned by the model and reconstructed whenever the davof scheme is
  (its fields stay unregistered; default check).
- `exact`: the control, the exact normal of the implicit function of the
  `exactNormal` sub-dictionary translated by U0 t at the centroid of the
  plane positioned by alpha along m_c (the cell centre when that plane does
  not exist); the state-time logic as before (default all).
- `ringArea` (new): n = sum_k m_k / |sum_k m_k| over the point-neighbour
  ring of the cell (the cells sharing at least one point with it, the cell
  included; zoneDistribute's stencil, across processor boundaries), only
  the interface cells contributing. This is the DAVOF normal of the union
  of the ring's cells: exact for a plane on any mesh, no division by a
  small area (default all).
- `ringFit` (new): a weighted least-squares quadric through the ring's
  recovered planes: in the frame of the cell's raw normal with the origin
  at the cell's own plane centroid, scaled by h_c, the Monge patch
  z = a + b x + c y + d x^2 + e xy + f y^2 fitted to the ring cells' plane
  centroids (one position equation each) and unit normals (two gradient
  equations each, as leia's `quadricFit` with `normalsOnly`), weights
  w_k = |m_k|/sum_j |m_j| (`weight area`) or alpha_k (1 - alpha_k)
  normalised (`weight alpha`); the normal of the fitted surface at the
  origin, (-b, -c, 1) rotated back. Ring cells with n_k.e3 <= `alignTol`
  (0.2) are dropped; fewer than `minCells` (6) usable ring cells or a
  singular system fall back to ringArea for that cell, counted as fit
  fallbacks (default all).

The deprecated keys `sliverNormal` (none | Youngs | plicRDF) and
`normalSource` (davof | exact | plicRDF) are mapped onto the models when
`normalModel` is absent, with a one-line Info note: sliverNormal Youngs ->
youngs/check and sliverNormal plicRDF -> plicRDF/check, both with
maxAngle 180 (the slivers-only rule they had), normalSource exact ->
exact/all, normalSource plicRDF -> plicRDF/all. Gate of the refactor:
the plane rows davof-detrixheAslam and davof-detrixheAslam-exact at
N = 20, Co 0.5, both arms, are bit-identical with the build before it
(`areaFractionError.dat`, `error.dat`, `volumeFractionError.dat` and the
written fields), through the deprecated keys and through `normalModel
cellLocal`; the plicRDF rows reproduce the cluster CSV to every digit
(max L1(n) 0.1420372339066551 at N = 20, normal arm).

Parallel check (2026-10-13): the plane case (normal arm, N = 20, Co 0.2)
with `normalModel ringArea`, `ringFit` and `cellLocal` on two processors
(`system/decomposeParDict`, simple (2 1 1), the cut through the interface
cells; `decomposePar`, `mpirun -np 2 advectorVoF -parallel`) gives the
same interface-cell counts per step as the serial run (193, 179, ..., 89)
and L1(alphaf) at round-off (1.7e-15 to 2.9e-15 against 1.8e-15 to
2.7e-15 serial), the normal error 0 to 2e-16: the ring stencil crosses
the processor boundary correctly.

Laptop measurements of 2026-10-13 (serial, Co 0.2, the new default of the
benchmark generators; CSVs `*_laptopCo0.2.csv` next to the cluster CSVs
of Co 0.5; the content bound and the second-layer bounding on, the
positionPlane fix in):

The plane (`planeTranslation3D`, ten steps, max over the steps of
L1(alphaf), the geometric-mean growth per step, the gate at 1e-12):

| row | arm | N | max L1(alphaf) | max L1(n) | growth (geo. mean) | gate |
|---|---|---|---|---|---|---|
| davof-detrixheAslam | normal | 20 | 3.37e-15 | 5.01e-14 | 1.07 | PASS |
| davof-detrixheAslam | normal | 40 | 5.70e-15 | 3.59e-13 | 1.07 | PASS |
| davof-detrixheAslam | oblique | 20 | 1.27e-14 | 6.41e-13 | 1.23 | PASS |
| davof-detrixheAslam | oblique | 40 | 1.56e-14 | 4.70e-13 | 1.10 | PASS |
| davof-detrixheAslam-exact | normal | 20 | 2.62e-15 | 1.97e-15 | 1.05 | PASS |
| davof-detrixheAslam-exact | normal | 40 | 3.68e-15 | 3.76e-15 | 1.00 | PASS |
| davof-detrixheAslam-exact | oblique | 20 | 8.08e-15 | 1.75e-15 | 1.06 | PASS |
| davof-detrixheAslam-exact | oblique | 40 | 7.21e-15 | 3.63e-15 | 0.94 | PASS |
| ringArea | normal | 20 | 2.72e-15 | 3.14e-15 | 1.05 | PASS |
| ringArea | normal | 40 | 3.70e-15 | 4.77e-15 | 1.00 | PASS |
| ringArea | oblique | 20 | 8.20e-15 | 1.02e-14 | 1.05 | PASS |
| ringArea | oblique | 40 | 7.94e-15 | 1.11e-14 | 0.95 | PASS |
| ringArea-check | normal | 20 | 3.43e-15 | 1.99e-14 | 1.07 | PASS |
| ringArea-check | normal | 40 | 5.12e-15 | 3.03e-14 | 1.05 | PASS |
| ringArea-check | oblique | 20 | 8.54e-15 | 3.26e-14 | 1.11 | PASS |
| ringArea-check | oblique | 40 | 8.80e-15 | 4.64e-14 | 1.01 | PASS |
| ringFit | normal | 20 | 2.70e-15 | 3.80e-15 | 1.05 | PASS |
| ringFit | normal | 40 | 3.55e-15 | 5.08e-15 | 0.99 | PASS |
| ringFit | oblique | 20 | 7.86e-15 | 1.43e-14 | 1.06 | PASS |
| ringFit | oblique | 40 | 7.51e-15 | 1.05e-14 | 0.94 | PASS |
| ringFit-check | normal | 20 | 3.32e-15 | 1.97e-14 | 1.06 | PASS |
| ringFit-check | normal | 40 | 4.92e-15 | 2.97e-14 | 1.05 | PASS |
| ringFit-check | oblique | 20 | 9.03e-15 | 3.33e-14 | 1.11 | PASS |
| ringFit-check | oblique | 40 | 9.35e-15 | 4.96e-14 | 1.01 | PASS |

Every row passes at Co 0.2, the own normal (cellLocal) included, whose
normal arm failed at Co 0.5 on the cluster (growth 6.8 and 23 per step at
N = 20 and 40) and still fails at Co 0.5 with this build (bit-identical:
1.36e-7 after ten steps at N = 20): the sliver amplification of the
normal arm is a Courant-number effect, absent at Co 0.2 over ten steps.
The ring rows replace every interface cell (208 at N = 20) and ringFit
falls back in 2 cells (the corner cells with fewer than six ring planes).

The sphere (`sphereTranslation3D`, N = 32, R = 4.8 h, 160 steps to t = 1;
L1(alphaf) after steps 1, 2, 3, 5, 10 and at the end; the end values of
the normal error L1(n), the centre error L1(x) and E1; min(alpha) over
the run):

| row | L1(alphaf) 1 / 2 / 3 / 5 / 10 / end | L1(n) | L1(x) | E1 | min alpha | cells (end) |
|---|---|---|---|---|---|---|
| davof-quadraticFaces | 4.39e-3 / 5.50e-3 / 8.06e-3 / 1.60e-2 / 6.15e-2 / 2.08e-1 | 1.29 | 1.54e-1 | 5.93e-3 | -9.79e-9 | 4910 |
| davof-quadraticFaces-exact | 4.38e-3 / 7.65e-3 / 5.71e-3 / 5.99e-3 / 6.65e-3 / 1.04e-2 | 3.64e-3 | 9.83e-4 | 2.76e-4 | -8.67e-15 | 426 |
| ringArea | 5.17e-3 / 6.12e-3 / 5.97e-3 / 6.98e-3 / 7.68e-3 / 1.78e-2 | 4.06e-2 | 1.35e-3 | 3.85e-4 | -3.40e-21 | 424 |
| ringArea-check | 4.43e-3 / 5.51e-3 / 6.28e-3 / 8.65e-3 / 1.15e-2 / 4.49e-2 | 1.32e-1 | 3.36e-3 | 8.13e-4 | -8.10e-15 | 487 |
| ringFit | 4.38e-3 / 7.71e-3 / 5.77e-3 / 6.04e-3 / 6.69e-3 / 7.93e-2 | 2.05e-1 | 1.32e-2 | 1.15e-3 | -2.45e-13 | 687 |
| ringFit-check | 4.39e-3 / 7.89e-3 / 8.48e-3 / 8.95e-3 / 1.17e-2 / 9.90e-2 | 3.01e-1 | 1.69e-2 | 1.23e-3 | -2.03e-9 | 765 |
| plicRDF | - | 3.92e-2 | 2.05e-3 | 3.18e-4 | -1.65e-24 | - |
| isoAlpha | - | 5.94e-2 | 2.17e-3 | 3.55e-4 | -3.31e-23 | - |

With the phantom-plane fix every row stays bounded. The own normal
(cellLocal) still loses the state on the sphere (L1(alphaf) 0.21, the
interface-cell count 416 to 4910, the normal error 1.29 at the end,
growth 1.2 per step over the first 20), as the plane gate at Co 0.5
predicted for a curved interface. The ringArea normal (applyTo all)
closes the loop: the interface-cell count stays at 416 to 428 over the
160 steps, L1(alphaf) grows from 5.2e-3 to 1.8e-2 (the exact-normal
control: 4.4e-3 to 1.0e-2), the normal error 4.1e-2 against plicRDF's
3.9e-2 and isoAlpha's 5.9e-2, the centre error 1.35e-3 against 2.05e-3,
E1 3.85e-4 against 3.18e-4, mass conserved to 3.5e-17. ringArea in check
mode (slivers and 15-degree outliers only, 48 to 167 cells replaced per
step) is three times worse in L1(n) and drifts; ringFit (all) matches
ringArea for 80 steps (1.35e-2 at step 80) and then loses cells (176
fit fallbacks and 672 interface cells at the end; the quadric fitted to
the drifting planes of the ring is the less robust of the two); the
check variants of both are worse than their all variants. The content
bound was active on 1845 faces (ringArea) to 331700 (cellLocal) over the
run. The N = 64 and 128 rungs and the orders are for the cluster.

The ring fallback of the plicRDF-normal control (2026-10-14, laptop,
sphere, Co 0.2; cases under `_laptop/sphere_fallback`). The control
`davof-<method>-plicRDF` (leia's `davof-plicRDFnormal`) does not converge
on the translated sphere in leia's hex study although plicRDF alone does
on the same alpha, and its Info line shows "kept m" cells (8 to 22 of
about 420 at N = 32, up to 42 of 1800 at N = 64, growing along the run):
the cells the davof scheme marks by their face fractions and the helper
plicRDF does not. The row `davof-plicRDFnormal-ring` (plicRDF, all,
`noNormalFallback ringArea`) gives those cells the ring normal. RESULT:
the row reproduces the keep control TO EVERY DIGIT, at N = 32 (L1(alphaf)
2.069546793547e-02, L1(n) 4.032040426541e-02, L1(x) 1.581827656378e-03,
Eshape 4.308063062444e-04, min alpha -3.83e-22; the keep control re-run
the same day gives the same digits and matches leia's committed hex
study) and at N = 64 (1.166392329698e-02 / 3.319389543382e-02 /
4.684655316210e-04 / 1.154731269617e-04, min alpha -9.99e-15, the
committed cluster values of the keep control), with the fallback taken in
exactly the former "kept m" cells (N = 32: 8 ring fallbacks at step 2, 22
at most; N = 64: 8 / 28 / 34 at steps 2 / 100 / 320, 42 at most). The
mechanism: the helper's criterion is alpha strictly inside (surfCellTol,
1 - surfCellTol), so a cell without a helper normal has alpha at 0 or 1
within 1e-8 while carrying face fractions above wispTol; `positionPlane()`
calls `vofCutCell`, which returns -1 or 1 for such an alpha before cutting
anything, so the cell gets no plane whatever normal it is handed. The
normal of a "kept m" cell never enters the transport; the fallback is
inert by construction, and the control's non-convergence has another
cause (open). Gate of the keyword: `ringArea` at N = 32 with the default
`keep` is bit-identical to the committed `_laptop/sphere_Co0.2` case
(`areaFractionError.dat`, `error.dat`, `volumeFractionError.dat`).

The benchmark generators compute the time step from `--Co` (default 0.2,
Scheufler and Roenby 2019 demonstrate second order for CFL 0.2 and
below; the cluster CSVs of 2026-10-13 are at 0.5), written into the case
(`CO`) and the results CSV (`Co`). Rows (`advectionRows.py`): the old
names map onto the models (`davof-<method>` cellLocal, `-Youngs<c>`
youngs/check/maxAngle 180, `-sliverPlicRDF<c>`, `-exact`, `-plicRDF`),
new `ringArea`, `ringArea-check`, `ringFit`, `ringFit-check` and, on
request (`--rows`), the control `davof-plicRDFnormal-ring` on the first
state construction given; every row carries `NONORMALFALLBACK` (default
`keep`) in `system/simulationParameter`, rendered as `noNormalFallback`
in `fvSolution` and reported in the results CSV. `getResults.py --cases
DIR --suffix S` collects another case tree into sibling CSV names.

### The mixed normal and the normal-source controls (2026-10-13, superseded by the normal models above)

The closed loop (own normal from alphaf, the alpha-positioned plane, alphaf
from the moved plane) amplifies through sliver interface cells (small
|m_c| compared with h^2: the positioned plane's sensitivity h/A, the own
normal's 1/A). Following Basilisk's mixed Youngs-centred scheme (Aulisa et
al. 2007, `src/myc.h`) the `davof` scheme keeps the slivers as interface
cells and gives them a stencil normal; keys of the alpha solver dictionary
(`davof.H`; since the normal models these keys are deprecated and mapped,
see above):

- `sliverNormal` (`none` | `Youngs` | `plicRDF`, default `none`) and
  `sliverTol` (c, default 0.1): a cell is a sliver when
  |m_c| < c V_c^(2/3). For the slivers the plane is positioned by alpha
  along the fallback normal: `Youngs` = `fvc::grad(alpha)` (points into
  alpha = 1 like m_c and `normal_`; gradAlpha's convention), `plicRDF` =
  the normal of a helper plicRDF scheme owned by the davof scheme,
  reconstructed on the same alpha whenever the davof scheme is (its
  `interfaceNormal`/`interfaceCentre` fields and its `reconstructionScheme`
  dictionary are refused by the registry, which keeps the davof scheme's).
  A sliver whose fallback normal is zero keeps m_c and is counted ("kept
  m"). The non-slivers keep the DAVOF normal. The Info line prints the
  sliver count per reconstruction; `sliverFlag.<group>` is written (1 in
  slivers).
- `normalSource` (`davof` | `exact` | `plicRDF`, default `davof`): the
  controls. `exact`: EVERY interface cell takes the exact normal of the
  implicit function of the `exactNormal` sub-dictionary (`functionType`,
  `origin` at t = 0, `U0`, the function's keys, as in `areaFractionError`)
  translated by U0 t, i.e. grad value at the centroid of the plane
  positioned by alpha along m_c (the cell centre when that plane does not
  exist), alphaf still transported by the moved plane; `plicRDF`: every
  normal from the helper plicRDF. The state time of `exact`: the first
  reconstruction of a time index is the advection step's and describes the
  state at t - dt, every later one in the same index (the function objects
  at the end of the step) the state at t, the constructor's the state at
  the start time; the time used is printed. The explicit position and the
  two plane defects are always formed from m_c.

Sanity on the plane (normal velocity, N = 20, ten steps, 2026-10-13;
`run/benchmark/advection/planeTranslation3D`): `exact` stays at round-off
for all ten steps (max L1(alphaf) 7.5e-15, max L1(n) 4.7e-15: the face
update is sound, the growth is the own-normal feedback); `davof` reproduces
the earlier history (1.4e-7 after ten steps); `Youngs 0.1` marks 38 of 208
cells and starts at 2.9e-4 (the Youngs normal of a plane is not exact) and
grows by about 2 per step to 0.19; the `plicRDF` normal starts at 4.7e-5
and ends at 2.2e-2. The sphere ladder on the cluster decides
(`run/benchmark/advection/run-benchmark.sbatch`).

The function object `areaFractionError` (`src/postProcessing`) measures the
transported face fractions against the exact DAVOF state of an implicit
function translated by `U0 t` (davofInit: `exactSphere`, `quadraticFaces`
or `detrixheAslam`), over the faces of the interface cells of either state
(L1, L2, Linf as fractions, the counts, the consistency
max |m_c - A n|/V^(2/3) of the transported state with the plane the scheme
positions in it), and writes the reconstruction error (normal, centre)
against the translated function at the end time (or every step). The
exact-solution gate is `run/benchmark/advection/planeTranslation3D`
(README there: the oblique arm passes at round-off; the normal arm is
exact per step and amplifies round-off through sliver interface cells by
a factor that falls with wispTol, a C2 finding).

### The cost of the DAVOF rows (2026-10-15)

The davof rows of the sphere study cost 1.16 times plicRDF per run
(leia `davof_advection_sphereTranslation3D_hex.csv`, N = 128, Co 0.2:
3797 s ringArea against 3283 s plicRDF) although the method's own work
(one closed-form sum per interface cell, one positioning) is below
plicRDF's iteration. The overhead sat in full-mesh passes repeated at
every reconstruction and every step, and in the diagnostic function
object. Removed without changing a floating-point result (every change
is argued for bit identity in the commit message; the cluster repetition
of the hex study against the committed tables is the check):

1. `davof::reconstruct()` positions each interface cell's plane ONCE,
   along the selected normal (above, "Normal models"); the exact and
   ringFit models keep the two passes because they read the planes
   positioned along m_c.
2. `davofInit::areaNormal` skips the faces with alpha_f = 0 (0 S_f adds
   +-0, which changes no sum, the sums start at +0 and +0 +- (+-0) = +0);
   m_c is formed in `areaNormalDavof` itself instead of a copy; the
   all-dry cells (m exactly +0, never touched) skip the all-equal test;
   the face moment of the explicit position is summed over the faces of
   the interface cells only, in ascending face order (the order of the
   former pass over all faces; `primitiveMesh::cells()` lists the owned
   faces first, so the cell's list is sorted first); the face-value mask
   is built once per mesh; the planes and the per-cell diagnostics are
   reset at the previous interface cells only (zero elsewhere by
   construction) instead of over the mesh.
3. The ring models: the exchanged field (m_c in the interface cells) is
   zeroed at the cells of the previous fill instead of over the mesh.
   The stencil and exchange set-up (`zoneDistribute::setUpCommforZone`)
   stays per call: `zoneCPCStencil` keeps every stencil it has computed
   (`uptodate_`), so a call costs one scan of the zone, the boundary
   points and the NEW interface cells only, plus in parallel one
   all-to-all for the pattern, which must follow the zone (it changes
   every step); caching it across the calls of one time index is unsafe
   because `zoneDistribute` is one object per mesh, shared with the
   helper plicRDF of the plicRDF model, whose own set-up with its own
   zone can run in between (the ring fallback of `davof-plicRDFnormal-
   ring`), and the gain would be one scan per extra reconstruction,
   which the `consistency` switch below removes anyway.

## The advection benchmarks (run/benchmark/advection)

- `planeTranslation3D`: the exact-solution gate of the face update (above).
- `sphereTranslation3D`: a sphere (R = 0.15) carried by U = (1 0 0) over
  one unit (3.3 diameters) in the 2 x 1 x 1 box at Courant 0.2 (the
  default of `genCases.py --Co`; the cluster CSVs of 2026-10-13 are at
  0.5), N = 32, 64, 128 cells across the height; rows isoAdvection + davof
  on the `initDavofFields` state (quadraticFaces, detrixheAslam) with the
  normal-model rows, + plicRDF and + isoAlpha on the setAlphaField state;
  `areaFractionError` every step, the reconstruction error at the end,
  `volumeFractionError` every step. README there for the measured
  numbers. With the phantom-plane fix every row stays bounded; the own
  normal (cellLocal) is bounded but loses the state on the sphere (the
  interface-cell count 416 to 4910 over the 160 steps at N = 32, the
  face-fraction error 0.21 at the end), while `ringArea` closes the loop
  (see "Normal models" above for the numbers; the earlier unbounded
  history of 2026-10-13, min(alpha) at -Co within a few steps, was the
  phantom-plane defect). The deformation regression
  (`testsuite/advection/test-deformationSphere` at 32 x 32 x 64, davof
  quadraticFaces, serial) predates the fix as well (Ebound 1.4 at t = 3,
  Emass 4e-3, Eshape 3.8e-2; plicRDF: 2e-8, 3e-16, 4.5e-3) and is to be
  repeated.

### Abort guard

`advectorVoF` stops a diverging run on request (keys of the top-level
`isoAdvector` dictionary of `system/fvSolution`, the one `period` and
`reverseTime` are read from; all off by default, so the testsuite is
unchanged):

```
isoAdvector
{
    abortEnabled              true;   // default false
    abortMinAlpha             -0.5;   // abort when min(alpha) < abortMinAlpha
    abortMaxAlpha             1.5;    // abort when max(alpha) > abortMaxAlpha
    abortInterfaceCellFactor  10;     // abort when the interface-cell count of
                                      // the reconstruction scheme exceeds this
                                      // factor times its count at the first
                                      // step; 0 disables this criterion
}
```

After every step's `t = ..., sum(alpha*V) ...` line the count of
`surf().interfaceLabels()` is reduced over the processors and the three
criteria are checked in that order. On a violation the solver prints one
line,

```
advectorVoF: ABORT at step <timeIndex>, t = <t>: <reason> (min(alpha) = <aMin>, max(alpha) = <aMax>, interface cells <n> of initial <n0>)
```

with `<reason>` one of `min(alpha) below abortMinAlpha`, `max(alpha) above
abortMaxAlpha`, `interface cells above abortInterfaceCellFactor times the
initial count`, writes the fields (`writeNow`) and ends with the normal
`End` and exit code 0 on purpose: leia's log classifier then reads the run
as COMPLETED and the workflow's converter records `diverged` from the
ABORT line. Measured on the sphere at N = 32, Co 0.2, row
davof-quadraticFaces (the own normal, 416 interface cells at the first
step; laptop, 2026-10-13): the factor 10 fires at step 111 (t = 0.69375,
4171 interface cells after 4141 at step 110; min(alpha) -6.25e-10,
max(alpha) 1), the run ends with "End" and the fields of t = 0.69375
written; the plane case without the keys (N = 20, normal arm,
davof-detrixheAslam-exact) is unchanged (max L1(alphaf) 2.61e-15 over
the ten steps, no ABORT line).

## How leia finds it

leia's `etc/leia-env.sh` appends `$LEIA_TPF_DIR/platforms/<WM_OPTIONS>/{bin,lib}` after
leia's own directories; the default is the `TwoPhaseFlow` clone next to the leia clone.
The `interFlow` benchmark arm of leia's studies (`config/interFlowDroplet{2D,3D}.yaml`)
therefore runs this build without any study-config reference to this repository.
