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

1. After the bounding every face no plane has set holds the upwind cell's
   new alpha, boundary faces the owner's (`defaultAreaFractions`; since
   2026-10-15 evaluated on the faces whose stored value can differ from
   that, see "The cost of the DAVOF rows", item 7). A sweep that cannot
   set its face (rule 7) writes the upwind cell's alpha at the start of
   the step, the value the former pass before the sweep gave every face.
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
7. A sweep that cannot set the face (`subFaceAreaEnd() < 0`) writes the
   default (the upwind cell's alpha at the start of the step), leaves
   the face unset for the second ring and the re-evaluation of rule 1,
   and is counted as a fallback; values outside [0, 1] by more than
   1e-12 are counted as clips and clamped.
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

### Face clipping models (2026-10-15)

What the face update stores on a face a plane has set is a
run-time-selectable model (`src/VoF/davof/faceClipping`, base
`davofFaceClipping`, keyword `faceClipping` of the alpha solver
dictionary, default `plane`; constructed by isoAdvection when the
reconstruction scheme is `davof`). `isoAdvection::storeAreaFraction`
hands the model, for every face a plane sets (the downwind faces of a
surface cell, the upwind faces its plane sweeps, the boundary faces, the
second ring; the processor sync exchanges the final values), the plane's
cell, the unit normal INTO the liquid, the plane's centroid translated to
the end of the step, the chord (the plane's trace on the face at the end
of the step, `cutFaceAdvect::chordEnd()`) and the plane-clipped area
(`subFaceAreaEnd()`), and stores what the model returns as the fraction.

- `plane`: the plane-clipped area itself, the face update as it was. Its
  `corrects()` is false, so isoAdvection keeps no chords and never calls
  the model: every existing row is bit-identical with the build before
  the models, by construction. The only additions on the path of a
  plane-set face are one tested bool (`clipCorrects_`, false), the split
  of `Un0 = UInterp.interpolate(x0, celli) & n0` into the vector `U0`
  and the same dot product `U0 & n0`, and the appends of `U0` to the
  per-plane lists; `cutFaceAdvect` fills `chordEnd_` only after
  `keepChord(true)`; the scheme's `usedNormal_` is a separate array; the
  tensor fields are allocated on request only. No floating-point
  operation of the sweep, the bounding or the reconstruction changes and
  no counter changes. The check on the cluster is the plane gate
  (`planeTranslation3D_hexCheck`-style, every value to all digits) and
  the hex sphere study's end values.
- `planeSagitta`: the sagitta correction, below.

The geometry and the sign convention (verified in the code's
conventions): the davof scheme's `normal_` and `m_c` point INTO the
liquid (alpha = 1; cutCell's `(p - base).n > 0` is liquid); isoAdvector's
plane normal `n0 = -normal/|normal|` points OUT of the liquid (its
`Un0 > 0` means "the cell fills"), so `storeAreaFraction` passes `-n0` to
the model. The curvature tensor `B` is the Hessian of the Monge patch
`z = z_0 + (1/2) x.B.x` of the interface over the positioned plane, `z`
along the normal `n` into the liquid, `x` in the plane from the plane's
centroid, in the frame of `davofNormalModel::mongeFit` (the ringFit
model's fit, factored out and shared; `e3 = n`, and for the ringFit model
`e3 = m_c/|m_c|` as before): for a liquid sphere the ring's plane
centroids lie on the liquid side of the tangent plane, `z = (x^2 + y^2)
/(2R)`, so the fitted `d, f` are positive and `B = I_t/R` is positive
definite, `tr B = 2/R`: "the surface bends toward the liquid" is
`t.B.t > 0`. `B` is stored as a symmetric tensor in the global frame
(`curvatureTensorDavof.<group>`, `t.B.t` for any in-plane `t`), with the
offset `z_0 = -(1/2) B : M` (`curvatureOffsetDavof.<group>`, `M` the cut
polygon's second moments about its centroid): the plane is positioned by
the cell's volume, so it is a SECANT of the surface over the cell and the
paraboloid's vertex sits at `z_0` (below the plane, `z_0 < 0`, for a
convex liquid region). Both fields are zero outside the interface cells
and in the fallback cells (fewer than `tensorMinCells` (6) usable ring
planes, `n_k.e3 <= tensorAlignTol` (0.2), a singular system; counted
"tensor fallbacks"); they are allocated on the first update, written with
`writeCurvatureTensor true` (default false). The tensor of the
start-of-step reconstruction serves the whole step (`advect` fits it right
after `reconstruct`; a translation does not change the curvature). The
scheme also exposes `usedNormals()`, the unit normal every interface
cell's plane was positioned along.

The correction. On a face with unit normal `n_f`, chord end points `p1,
p2` (length `L`, direction `t`), the surface's trace is displaced from the
chord in the face plane, along the chord's in-plane normal into the
liquid, by `w/cos(phi)`, where `w(l) = z_0 + (1/2) x(l).B.x(l)` is the
surface's offset from the plane along `n` at the chord point `x(l) = a +
l t`, `a = p1 - c'`, `c' = x0 + U0 dt` the centroid at the end of the
step, and `cos(phi) = |n - (n.n_f) n_f|` (Meusnier: the trace's curvature
is the normal curvature over `cos(phi)`). Where `w > 0` the surface lies
on the liquid side of the plane and the true liquid area is smaller than
the plane's:

    dA = -(1/cos(phi)) int_0^L w(l) dl
       = -(1/cos(phi)) [ z_0 L + (1/2)( a.B.a L + a.B.t L^2 + t.B.t L^3/3 ) ]

exact to O(h^4) in the area for a second-order normal and a first-order
`B` (the end triangles where chord and arc meet the face's edges are
O(h^4); an O(h) error of `B` moves the area by O(h^4)), so the stored
fraction is second order where the plane clipping is first order
(`sagittaForm full`, the default). The original formulation, the
parabolic segment between the chord and the arc through its end points,
`dA = +(2/3) L s / cos(phi)` with the sagitta `s = t.B.t L^2/8`, i.e.
`t.B.t L^3/(12 cos(phi))`, is kept as `sagittaForm segment`: it assumes the
chord's end points ON the arc, which the volume-positioned plane does not
give, and drops the offset and the `a`-terms, which are of the same
order. The check on a cube of side `h` cut through its centre by a plane
normal to `e3` with `B = diag(B11, B22)`: `M = diag(h^2/12, h^2/12)`,
`z_0 = -(B11 + B22) h^2/24`; on the face `x = h/2` (`L = h`, `t = e2`,
`a = (h/2, -h/2)`) the full form gives `int w = B11 h^3/12`, `dA = -B11
h^3/12`: the plane OVERESTIMATES the wetted area of every side face of a
convex liquid region (the surface's trace lies on the liquid side of the
secant's trace), while the segment form gives `+B22 h^3/12`, the opposite
sign and the curvature of the other direction. Guards (keywords, with
their defaults): no correction without a chord (the face on one side of
the plane at the end of the step: a corner triangle of O(h^4) at most)
or in a tensor-fallback cell; `sagittaMinCos 0.2` (`cos(phi)` below it:
the interface nearly parallel to the face, the trace ill-conditioned;
"skipped (cos)"); `sagittaMinChord 0.05` (`L` below it times
`sqrt(|S_f|)`: the plane grazes a corner; "skipped (chord)");
`sagittaMaxFraction 0.25` (`|dA|` capped at it times `|S_f|`; "capped");
the stored area clamped to `[0, |S_f|]` ("clamped"); `tensorSource ring`
(the ring Monge fit, the only source in this round). Planar faces: the
chord is the plane's trace at `t = dt` on the face's edges (the cut
points of the last quadrilateral of `timeIntegratedArea`, the cut points
of `calcSubFace` for a stationary plane); on a non-planar face (the fan)
the crossings of the plane with the face polygon's edges and the face's
mean normal are used, an approximation (the perturbed-mesh transport is
a separate matter). Configuration (the alpha solver dictionary; the old
configuration without the keyword is the default):

```
faceClipping        planeSagitta;   // plane (default) | planeSagitta
sagittaMinCos       0.2;
sagittaMinChord     0.05;
sagittaMaxFraction  0.25;
sagittaForm         full;           // full | segment
tensorSource        ring;
tensorAlignTol      0.2;            // the scheme's fit (davof keys)
tensorMinCells      6;
writeCurvatureTensor false;
```

Info lines: the scheme prints `davof: curvature tensor on N cells, F
tensor fallbacks` after the start-of-step reconstruction, and the face
update's line gains `; sagitta corrected F faces, C capped, K clamped, S
skipped (cos), M skipped (chord)` (only when the model is not `plane`;
leia's `tpf_advection_metrics.py` parses both).

Expected on the plane test (`planeTranslation3D`): the ring's planes of
an exact plane state have centroids on the plane and exact normals to
round-off, so the fitted quadratic coefficients are round-off (`d, e, f`
about 1e-15 in the scaled units, `B` about 1e-15/h, 5e-12 1/m at N = 20
of the 4 mm box) and `|dA|/|S_f|` about `B h/12` ~ 1e-16 per face, the
offset term `z_0 L/|S_f|` ~ `B h/24` as small: the plane gate at 1e-12
must hold for the sagitta rows as for the plane rows (a tensor of the
order of 1e-13 1/m moves a fraction by 1e-13 h/12, far below the gate at
every N). Expected on the sphere: the face-fraction error second order
with the exact normal (the clipping law's q = 1 was the plane's, the
corrected fraction is O(h^2)), the own normal recovered at second order
where the slivers take the ring normal; the ring normal's first-order
bias on a curved interface limits its row (the slope error O(h) enters
`w` at O(h^2)). The sphere's `z_0`: the plane positioned by alpha in a
cube cut through its centre, `z_0 = -h^2/(12 R)`, 1.7e-3 h at R = 4.8 h.

### The curved clipping surface and its volume rules (2026-10-23/24)

`clippingSurface exactSphere | paraboloid` in the isoAdvector dictionary replaces the plane
sweep by ONE quadric surface per interface cell for both the volume flux (the liquid volume of
the prism swept between a face's pre-image f - U dt and the face under the upwind cell's
surface) and the face update (the pre-image's liquid fraction under the same surface), so that
alpha_k and alpha_f are the moments of one surface at every step (`isoAdvectionCurved.C`,
`davofCurvedGeometry.H`). `exactSphere` takes the sphere of the `exactNormal` dictionary at the
state's time (the machinery gate for the sphere translation case); `paraboloid` the paraboloid
in the frame of the normal the scheme used, with `fitDof six | one | none` (how many of its
coefficients the cell's own face fractions set) and `fitAnchor faces | volume` (the offset from
the faces or from alpha_k V). The stable variant is `fitDof none; fitAnchor volume`: the frame
from the normal, the curvature tensor from the ring, the offset from the conserved alpha_k.
Serial runs only until the surfaces are exchanged across processor boundaries.

The three geometric quantities, the cell volume under a surface (the anchor, the surface
defect), the prism volume (the flux) and the face fraction (the update), come from one
run-time selected rule, `volumeRule` (`davof/volumeRules`, `davofVolumeRule`, RTS as the normal
models):

- `splitQuadrature` (the default, the reference): the triangle integrals split at the exact
  traces of the quadric (the roots along every sweep line and the swept edge, three evaluations
  each), 12-point Gauss-Legendre per piece, refined adaptively in both sweep directions to
  `curvedTolRel` of the volume (default 1e-11, depth cap `curvedMaxDepth`). Any quadric. On the
  translated plane at N = 40 it reproduces the exact linear cut to 2e-14 of the cell.
- `kromerBothe`: J. Kromer, D. Bothe, Third-order accurate initialization of volume fractions
  on unstructured meshes with arbitrary polyhedral cells, J. Comput. Phys. 475 (2023) 111840,
  doi:10.1016/j.jcp.2022.111840, arXiv:2111.01073 (building on J. Comput. Phys. 396 (2019)
  761-784, doi:10.1016/j.jcp.2019.07.005, arXiv:1805.03136). Their local paraboloid is the DAVOF
  surface itself: the divergence theorem with x - o gives the face sum plus the paraboloid
  patch integral, which Green's theorem turns into the line integral of H1 = c t1 - (b11 t1^3/6
  + b12 t1^2 t2/2 + b22 t1 t2^2/2) along the conic traces on the faces; the wetted face areas
  by Green's theorem as well. The conic is classified in principal form (their Eqs. 26-34) to
  order the boundary crossings and to split long arcs; every arc is integrated as the graph
  over its chord with the exact root along the chord's normal, adaptive Gauss-Legendre 12
  (`davofKromerBothe.H`). A planar face (every vertex within 1e-10 L_s of the face plane) is
  integrated whole, a warped face as the fan triangles about the face centre; the key
  `kromerBotheFan` overrides the test (`auto`, the default, the test; `true` every face
  fanned, the rule before 2026-10-25; `false` every face whole). The same choice applies to
  the cell faces of the volume anchor, the top and bottom of the flux prisms and the displaced
  faces of the face fractions. On the hexahedral sphere the whole faces cut the step by 45 %
  at N = 64 and 40 % at N = 128, the end error changing in the eighth digit (L1(alphaf)
  0.00195336 against 0.00195339). Falls back to `splitQuadrature`
  per polyhedron where a configuration is not resolved, and for a sphere surface; counted in
  the Info line ("volume rule kromerBothe fallbacks N") and, when any occur, by reason
  ("kromerBothe fallbacks by reason: conic, point, root, gradient, area, volume"). On a flat
  arc (a line trace, or a parabola of round-off curvature, as on the translated plane where
  the ring tensor is 1e-10) the arc's midpoint lies on its chord to round-off and the root
  along the chord's normal is the one nearest to the chord; selecting it by the arc's side
  took the far root of the near-degenerate quadratic in up to 70000 prisms per step and
  sent them to the capped fallback (fixed 2026-10-25; the plane test is at round-off and
  without fallbacks since).
- `exactSphere`: S. Strobl, A. Formella, T. Pöschel, Exact calculation of the overlap volume of
  spheres and mesh elements, J. Comput. Phys. 311 (2016) 158-172, doi:10.1016/j.jcp.2016.02.003,
  through the authors' reference implementation (`davof/thirdParty/overlap/overlap.hpp`, MIT,
  header-only on Eigen 3; `src/VoF/Make/options` adds `-isystem $(EIGEN_INC)`, which
  `etc/davof-env.sh` sets to `$HOME/OpenFOAM/thirdparty/eigen/eigen-3.4.0` when unset). The
  cell as the (centre, face centre, edge) tetrahedra of every face's fan, the prism as three
  tetrahedra per fan triangle; the face fraction by the quadratic face cut (exact for the
  circle trace up to the sliver's 12-point rule). For `clippingSurface exactSphere`; a
  paraboloid falls back to `splitQuadrature`, counted.

Other algorithms of the curved flux, with their sources: the volume anchor's root (key
`anchorMethod`: `ccs`, the default; `newton`; `illinois`; `compare`; 2026-10-25). The slope
dV/dc of the cell's liquid volume in the offset is the area of the wetted patch projected along
the frame's axis, which the Kromer-Bothe rule returns as the line integral oint t1 dt2 around
the trace arcs it walks for the volume anyway (`cellLiquidVolumeSlope`). `ccs`: the
Consecutive Cubic Spline positioning of T. Marić, Iterative volume-of-fluid interface
positioning in general polyhedrons with Consecutive Cubic Spline interpolation, J. Comput.
Phys.: X 11 (2021) 100093, doi:10.1016/j.jcpx.2021.100093, arXiv:2012.02227, applied to the
paraboloid: once the root is bracketed (a Newton step from the start), the iterate is the root
of the cubic Hermite interpolant through the two bracket ends with the exact volumes and
slopes, re-bracketed at every step, an iterate outside the bracket or a bracket not halved in
two steps replaced by the bisection point; between the vertex offsets (the breakpoints of
V(c), where its second derivative jumps) V is smooth, a cubic for a plane, so the model is
exact there up to the quadrature. `newton`: Newton on the slope, safeguarded by the bracket
the iterates build (an iterate outside it is replaced by the bisection point: the rtsafe of
Press, Teukolsky, Vetterling, Flannery, Numerical Recipes, 3rd ed., Cambridge University
Press, 2007, Sect. 4.6 and 9.4, ISBN 978-0-521-88068-8) and corrected by the quadratic model
whose second derivative is the secant of the two latest slopes (Halley's step with that
estimate: T. R. Scavo, J. B. Thoo, On the geometry of Halley's method, Amer. Math. Monthly 102
(1995) 417-426, doi:10.2307/2975033). Both converge when the bracket closes to 1e-11 h or the
Newton step r/A is below it with the residual within `anchorTolRel`; where the slope is not
available (a fallback polyhedron, another rule) and no bracket exists yet, or with
`anchorMethod illinois`, the root by a bracket of doubling steps and the Illinois rule,
M. Dowell, P. Jarratt, A modified regula falsi method for computing the root of an equation,
BIT 11 (1971) 168-174, doi:10.1007/BF01934364 (the only method until 2026-10-25: 6.2
cell-volume evaluations per interface cell per step, 57 % of the advection step at N = 64).
Where the target is 0 or V (a full or empty cell flagged as an interface cell) V(c) is flat on
an interval of roots; Illinois stops at the first bracket end of its doubling steps (0.75 h or
1.75 h from the start) and the ladder's ring row carries that choice; the offset of such a cell
is DEFINED (2026-10-25) as the edge of the interval, the surface tangent to the cell (the first
offset at which the cell is entirely liquid or entirely empty), computed exactly as the extremum
of the offset function over the cell's vertices, edges and face planes (`kbOffsetRange`) and
moved 1e-12 h onto the flat side; `ccs` and `newton` take it (the Info line counts these cells,
"edge N"), `illinois` keeps its bracket end. `compare` runs Illinois and ccs from the same
start, continues with ccs, counts the cells whose offsets differ by more than 1e-9 h (and those
of them at the edge) in a second Info line and prints the first 20 per step. The Info line
counts the anchor's evaluations ("evaluations N (per cell x)"). The anchor's two face-fraction
diagnostics, the face residual before the anchor and the face shift after it (14 % of the
step), are computed and printed only with `anchorDiagnostics true` (default false). The
six-coefficient fit by
Gauss-Newton with the damping of D. W. Marquardt, An algorithm for least-squares estimation of
nonlinear parameters, J. Soc. Indust. Appl. Math. 11 (1963) 431-441, doi:10.1137/0111030; the
six-point degree-4 triangle rule of the retired unsplit quadrature from D. A. Dunavant, High
degree efficient symmetrical Gaussian quadrature rules for the triangle, Int. J. Numer. Meth.
Engng. 21 (1985) 1129-1148, doi:10.1002/nme.1620210612; the Gauss-Legendre nodes by Newton on
the Legendre polynomials (Press, Teukolsky, Vetterling, Flannery, Numerical Recipes, 3rd ed.,
Cambridge University Press, 2007, Sect. 4.6, ISBN 978-0-521-88068-8); the tetrahedral cell
fractions of the initialisation from M. Detrixhe, T. D. Aslam, From level set to volume of
fluid and back again at second-order accuracy, Int. J. Numer. Meth. Fluids 80 (2016) 231-255,
doi:10.1002/fld.4076. The closed forms of F. Evrard, R. Chiodi, A. Han, B. van Wachem, O.
Desjardins, First moments of a polyhedron clipped by a paraboloid, SIAM J. Sci. Comput. 45
(2023) A2250-A2274, doi:10.1137/22M1524308, arXiv:2210.07772, are the alternative to the
Kromer-Bothe rule for the paraboloid and are not implemented.

The volume flux of a face (key `fluxIntegration`, 2026-10-25): the phase-specific flux is phi_f
times the time integral of the wetted face fraction over the step. `prism` (the default): the
liquid volume of the prism between the face's pre-image f - d and f under the upwind surface,
that integral exactly for a uniform translation, by the volume rule. `trapezoidal`: phi_f dt
(alpha_f^n + alpha_f^{n+1})/2 from the face fractions the loop has, alpha_f^n the stored
fraction and alpha_f^{n+1} the pre-image's liquid fraction under the same surface, which the
face update computes in the same pass; no prism volume. `simpson`: phi_f dt (alpha_f^n + 4
alpha_f^{n+1/2} + alpha_f^{n+1})/6 with alpha_f^{n+1/2} the liquid fraction of the half-step
pre-image f - d/2 under the same surface (one more face fraction per face instead of the prism's
faces). The fraction of a face is a cubic in time between the passages of the trace over the
face's vertices, so Simpson's rule is exact on a step without such a passage and the
trapezoidal rule is second order. `simpsonSplit`: Simpson's rule on the sub-intervals between
those passages, whose times are predictable: the roots sigma in (0, 1) of phi(x_v - sigma d) = 0
along the trajectory of every vertex of the face over the step (and of the face centre, the apex
of a warped face's fan), a quadratic in sigma for the paraboloid; each sub-interval takes the
pre-image fractions at its ends and its midpoint (two more face fractions per sub-interval, and
few faces have a breakpoint in a step), so that on the plane the rule is exact in every step; the
time-quadrature analogue of the split quadrature's breakpoints. Where three or more vertices
meet the surface at one time the face lies on the surface then (a face parallel to a plane that
passes it: an axis-aligned plane on a hexahedral mesh; not the plane test, whose plane is
oblique to the mesh) and its fraction jumps from 0 to 1 or back; the end values of such a
breakpoint are the one-sided limits, evaluated 1e-9 of the step before and after it (the
pre-image is then entirely on one side). `gaussSplit`: the three-point Gauss-Legendre rule on
the same sub-intervals (exact for a quintic in time); its nodes lie inside the sub-intervals, so
no face fraction is evaluated at a breakpoint, the instant at which a vertex of the pre-image
lies exactly on the surface (the degenerate configuration of the face integral), and the stored
fraction alpha_f^n does not enter: every value is a fraction under the same surface, as the
prism volume is; three face fractions per sub-interval. The signs follow the prism path
(the flux carries the sign of phi_f); a face that does not sweep (|d . S_f| below 1e-10 |d| |S_f|)
gets no flux in every mode. The leia token is `TPF_FLUX_INTEGRATION` (default `prism`).

The leia token is `TPF_VOLUME_RULE` (default `splitQuadrature`; `cases/davof/*Translation3D`).

### Warped faces: the wetted area vector as the face state (branch davof-warped, 2026-10-25)

On a non-planar face the wetted part's area vector a_f = sum_t A_wet,t n_t (the fan triangles
(x_f, p_i, p_{i+1}) about the face centre, each with its own unit normal n_t) is not parallel to
S_f, and alpha_f S_f does not carry the normal identity sum_f a_f^out = -A n. Every fan triangle
contains x_f, so (x - x_f).n_t = 0 on it and the first moment of the wet part is
int_f chi (x - x_c).n dS = (x_f - x_c).a_f: the vector alone (three scalars) is the transportable
face state; no separate first moment is needed for the explicit position. On a planar face
a_f = alpha_f S_f to round-off.

Keys (the alpha-solver dictionary; every default is the behaviour before the keys, and every new
branch sits behind a null pointer or a false flag, so the default runs are bit-identical):

```
faceState            scalar;      // | vector: wettedAreaDavof.<group> (surfaceVectorField, m^2)
                                  //   transported next to alphaf; needs a clippingSurface
tensorCentroid       cut;         // | anchored: the centroids of the ring Monge fit of the
                                  //   curvature tensor re-anchored by the volume rule
kromerBothePrismApex vertexMean;  // | faceCentre: the warped flux prism's top and bottom fanned
                                  //   about the mesh face centre (the cell's own fan)
```

- Initialisation: `initDavofFields` with `faceState vector` in `initAlphaFieldDict` writes
  `wettedAreaDavof.<phase>` from the same fan cuts as alphaf (detrixheAslam: the linear cut's wet
  polygon per triangle, `davof::triNegativePolygon`; quadraticFaces: the quadratic fraction times
  the triangle's area vector; exactSphere: alpha_f S_f, its faces are planar); `invert` maps it to
  S_f - a_f; the consistency against the tet sum is printed for the vector state too.
- Face update (`isoAdvectionCurved.C`): `davofVolumeRule::faceWettedArea` returns the
  pre-image's wetted area vector with the fraction (Kromer-Bothe: `kbFaceWettedArea`, the Green
  integral's wetted area per fan triangle, or of the whole face where it is planar, times its unit
  normal; the base rule: the quadratic triangle cut per triangle, `faceCutDisplacedVector`). The
  volume flux is unchanged. Faces no surface sets take alpha_up S_f (`setDefaultAreaFraction`).
  Processor patches: the vector travels with alphaf in `syncProcPatches`, with its sign flipped
  (each side's patch face is oriented outward from its own cell, as the flux).
- Reconstruction (`davof.C`): m_c = sum_f a_f^out over the same wet-face list in the same order,
  and the face moment of the explicit position (x_f - x_c).a_f. Everything downstream (the normal
  models, the ring normal, the tensor fit, the positioning by alpha) reads m_c as before.
- The tensor centroids. OpenFOAM's `cutCellPLIC` (the positioning by alpha) closes the wet part of
  a warped face by the chord between the two edge crossings, while the fan face's wet region ends
  at the plane's trace on the fan triangles: the cut is that of another polyhedron than the fan
  polyhedron of the mesh, of the Detrixhe-Aslam alpha and of the volume rule, so the positioned
  plane and its centroid are off the true plane by a fraction of h even for exact data, and the
  ring Monge fit through those centroids returns a spurious tensor (a curved surface for a plane).
  `tensorCentroid anchored` re-anchors, before the fit, every interface cell's plane in the frame
  of its used normal through the cut centroid by the volume rule's cell volume (bracket over the
  cell's extent along the normal, the Illinois rule of Dowell and Jarratt, BIT 11 (1971) 168-174,
  doi:10.1007/BF01934364, to a residual of 1e-15 V or an offset bracket of 1e-13 h) and hands
  o + c e3 to the fit (`davof::tensorCentres()`, exchanged across processors in place of
  `centre_`); failures are counted in the Info line ("anchored centroids: failed N"). The scheme's
  `centre_` (the surface origin, the exact model's evaluation point, the reconstruction error) is
  unchanged; the offset of the surface is the volume anchor's in any case.
- The prism apex. The Kromer-Bothe rule fans a warped prism's top and bottom about their vertex
  mean, the cell's faces and the face fractions about the mesh face centre; on a warped face the two
  fans differ by a thin lens, whose liquid enters the flux and not the cell's own volume.
  `kromerBothePrismApex faceCentre` uses the mesh face centre and its translate.
- The flux prism over a warped face (two fixes found on the perturbed plane's oblique arm,
  U = (1 0 0), where every row and face state failed at 2e-3 to 1e-2 while the normal arm was
  exact; unconditional, the faces of a prism over a planar face unchanged vertex for vertex
  (there the interior-point test and the sign of d.S_f agree) and the hexahedral gate
  bit-identical): (1) `kbPrismFaces` oriented the prism's fan triangles and sides by an
  interior-point test, which on a thin prism over a warped face (d nearly parallel to the face)
  reads the warp; they are now oriented by the sign of d.S_f (`kbOrientedFan`), and the
  split-quadrature fallback `faceSweptLiquidVolume` gives a fan triangle that d sweeps backwards
  its negative sign. (2) The net liquid volume through such a face lies in [-V_back, V_fwd] (the
  backward and forward triangle sweeps), not in [0, |d.S_f|]; the rule's range check and clamps
  use that range (`kbPrismRange`) where a triangle is swept backwards.
- Measured on the perturbed hexahedra (every vertex displaced by up to 0.1 h; leia
  `planeTranslation3D_perturbedWarped`, Co 0.2, N = 20/40/80, both arms, the volume-anchored
  ring-tensor paraboloid with the Kromer-Bothe rule, `kromerBothePrismApex faceCentre`): with
  `faceState vector` and `tensorCentroid anchored` the exact-normal row is at round-off on both
  arms at every N (run maxima of L1(alphaf) 1e-14 to 6e-14, L1(n) 1e-14), the ring-normal row at
  round-off on both arms at every N as well (before the davof-sliver merge f133781 the oblique arm
  at N = 80 had isolated single-face events, L1 up to 1.1e-8, at faces with a vertex 2e-4 h from
  the start-of-step plane: the dropped sliver arc that merge removes; on f133781 the case is at
  7.0e-14). The scalar state fails in the
  ring row (L1(n) 2.1e-3 to 2.4e-3, L1(alphaf) 2.6e-4 to 4.1e-4, flat in N); `tensorCentroid cut`
  fails in every row at 1.3e-5 to 1.7e-5, flat in N; `kromerBothePrismApex vertexMean` fails
  at 2e-5 to 3.5e-5 with everything else in place. On the hexahedra `faceState scalar` is
  bit-identical with the unmodified build (8a1c5e4, then f133781 and 1f914f9 after the rebases
  onto the sliver and the perf merges) in every number of areaFractionError.dat, error.dat and
  volumeFractionError.dat (8 cases), `vector` at round-off; on 1f914f9 the perturbed plane's
  vector + anchored cases are at 5e-15 to 6e-14 in L1(alphaf) and 2e-16 in L1(n) throughout.
- Diagnostics: `areaFractionError` takes m_c of both states from the vectors when
  `wettedAreaDavof` is registered (the interface-cell set of the norms) and logs the vector error
  |a_f - a_f,exact|/|S_f| (L1, Linf); the .dat columns are unchanged. Its `consistency` column
  compares m_c with the area vector of the CUT polygon, which carries the chord closure above on a
  warped mesh, with either state.

The leia tokens (both `cases/davof/*Translation3D`, defaults as above): `TPF_FACE_STATE` (also in
`initAlphaFieldDict`), `TPF_TENSOR_CENTROID`, `TPF_KB_PRISM_APEX`.

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
against the translated function at the end time (or every step). Its cost
is a diagnostic's (2026-10-15): the exact state is built over the whole
mesh at every write (exactSphere: every face against the ball, every cell
box), two area sums follow, and the consistency column costs one EXTRA
reconstruction of the scheme per write. The switch `consistency` (default
true) turns that reconstruction off (the column reads -1); everything is
computed in `write()` only (`execute()` does nothing), so `writeControl
onEnd` measures the end state at no cost during the run, while
`writeControl timeStep` with an interval the end step does not divide
writes NOTHING (no `end()` of its own). The skipped reconstruction leaves
the transport bit-identical for the models that reconstruct from alphaf
and alpha alone (cellLocal, youngs, exact, ringArea, ringFit: the next
step's reconstruction is the first of its time index whatever happened
before, and the exact model's state time follows the index, not the
count); it CHANGES the path of the plicRDF model (the helper starts its
iteration from the previous normals), as `reconstructionErrorEveryStep`
does for the plicRDF comparator. The
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
4. `isoAdvection`: the second-ring sweep collected its candidates in two
   arrays of the size of the mesh allocated and filled per step and
   visited every face in a final pass; the arrays are now sized once,
   the candidate faces are listed as they are first reached, visited in
   face order (the order of the former pass) and reset afterwards. A
   face reached again by the same plane (through a second second-ring
   cell holding it) is skipped: its distance recomputes to the value it
   holds and cannot beat itself under the strict comparison, so the
   closest-plane choice is unchanged (same surface-cell order, the first
   of equal distances wins). The per-face records of `planeFlux`
   (`planeFluxUpCell_`, `dVfUpwindPlaneFace_`) and the cell-to-plane map
   are sized once and reset at the entries the step set (a list kept by
   `planeFlux`, `surfCells_`); `boundCellContent` marks its cells from
   that list instead of scanning every face (the bitSet is visited in
   cell order whatever the insertion order, so the sequential dependence
   between neighbouring bounded cells is unchanged).
5. `areaFractionError`: the switch `consistency` and `writeControl onEnd`
   (above).
6. Profiling: OpenFOAM's markers (`addProfilingInFunction`, free when
   profiling is off: the trigger constructor returns at once) sit in
   `davof::reconstruct`, `davofNormalModel::select` and
   `ringAreaNormals`, every model's `modelNormals` (one name per class:
   the marker is the function's signature), `isoAdvection::advect`,
   `timeIntegratedFlux`, `defaultAreaFractions`, `sweepSecondRing`,
   `boundCellContent`, `limitFluxes`, `boundFlux`,
   `applyBruteForceBounding`, `davofInit::areaNormal` and the three
   state constructions, `areaFractionError::write` and
   `volumeFractionError::write`, besides plicRDF's own. Switch them on
   with `profiling { active true; }` in `system/controlDict` (Time reads
   it; `cpuInfo true; sysInfo true;` are optional); the tree of
   cumulative times per marker (calls, total and children time) is
   written at every field write to `<time>/uniform/profiling`, so the
   last write of the run holds the whole run. The leia timing study
   `config/davof/sphereTranslation3D_hexTiming.yaml` (token
   `TPF_PROFILING`) uses it.

The profile of that study (sphere, N = 128, Co 0.2, 640 steps, serial,
the ringArea row; totals over the run) after items 1-6: `advect` 431 s;
`davof::reconstruct` 94 s, of which `areaNormal` 32 s (against
plicRDF's `reconstruct` of 265 s); `timeIntegratedFlux` 61 s, of which
`defaultAreaFractions` 17.6 s and the second-ring sweep 7.4 s; a second
`defaultAreaFractions` call from `advect` 19.2 s; `limitFluxes` 66.6 s
(plicRDF: 57 s). The second round (2026-10-15, commits c1b3995, 660d741
and bfe0c00) removes the remaining full-mesh
passes of the face update, each argued for bit identity (the cluster
repetition of the hex study against the committed tables is the check):

7. The default face fractions. The pass before the sweep (every face
   the upwind cell's alpha) is gone: its value was read only through
   the fallback of `storeAreaFraction` ("keeps the default", rule 7)
   and, in parallel, through the sync of such a fallback face;
   `sweepSecondRing` reads the mask alone, `planeSpeaks`, `planeFlux`
   and `boundCellContent` read alpha and `dVf`, and the reconstruction
   runs before the step. The fallback now writes that default itself
   (alpha at the start of the step; `alpha1_` is updated after the
   sweep), stays unset and is listed. The pass after the bounding is
   restricted to the faces whose stored value can differ from the
   target (alpha of the upwind cell after the bounding). INVARIANT
   between the steps: an unset face holds `alphaDefault_` of its upwind
   cell (boundary faces: the owner), `alphaDefault_` being a copy of
   alpha at the last evaluation, as long as nothing has written the
   face since; and every write of a face in a step (a plane, the sync,
   a fallback) marks or lists it. Hence the target can differ only on
   (a) the faces whose upwind direction changed, re-scanned when `phi`'s
   event number moved since the last evaluation (every non-const access
   of a GeometricField, `ref()`, `primitiveFieldRef()`,
   `boundaryFieldRef()`, `operator=`, `correctBoundaryConditions()`,
   records one; the translation cases of `advectorVoF` write `phi` at
   most once, at a reversal, the periodic and spiralling ones every
   step), (b) the faces set in the
   previous step (listed), (c) the fallbacks of this step (listed), (d)
   the faces of the cells whose alpha differs BITWISE from
   `alphaDefault_` (one `memcmp` per cell; -0 and +0 are different
   values in the written field, so `!=` would not do). Evaluating a
   face always writes the target, so a face visited twice is unchanged.
   The first evaluation, and the one after a topology change or a size
   mismatch, is the full pass and records the copy and the directions.
   What remains per step: one compare per cell and the faces of the
   changed cells.
8. The area normal. `isoAdvection` keeps a bitSet of the faces whose
   fraction is not 0, updated at every write of a face value
   (`setDefaultAreaFraction`, `storeAreaFraction`, the proc sync;
   rebuilt from the field in the full evaluation), and hands its sorted
   list (`bitSet::toc()`, ascending) to the `davof` scheme at the end
   of every step (`publishWetFaces`, `davof::setWetFaces`).
   `davof::reconstruct` (`formAreaNormal`) sums m_c over that list when
   it is valid (present, built for the current face count, no changing
   topology): the former pass over every face in face order skipped the
   faces with alpha_f = 0 (item 2), and the list holds exactly the
   faces with alpha_f != 0 in the same order (the internal faces, then
   the patch faces in patch order), so every cell receives the same
   terms in the same order; a face at 1 contributes 1 S_f as before,
   and the all-equal rule still zeroes a cell whose faces are all 1.
   The cells the sum touched are kept in a bitSet: they are zeroed
   before the next sum (every other cell holds the +0 of the last full
   zeroing, never written since), and they are the cells the
   interface-cell selection visits, ascending as the former loop over
   all cells, which skipped the cells with m = 0 (the untouched ones
   hold exactly +0). Without a valid list (the constructor, the first
   advection step, `reconstructInterface`, a topology change) the full
   pass runs and records the touched cells. The per-face copy
   `alphafAll_` is gone; the all-equal test and the face moment read
   the face value through a patch lookup, the same double.
   `davofInit::areaNormal` (the full pass with the per-face copy) is
   unchanged: `maxConsistency` and `areaFractionError` use it for the
   exact state and for the transported field.
9. `extendMarkedCells` (the one-layer extension of `limitFluxes`, called
   once by upstream and a second time for the second layer, rule 9):
   the cells were marked by testing every face of the mesh against the
   marked-face bitSet; now the set bits of that bitSet are visited (the
   iterator skips the zero words). The resulting SET of cells is the
   same, and its consumers (`boundFlux`, the second-layer count) iterate
   the bitSet in cell order whatever the order the bits were set in.
   This is where the extension's 10 s sat (a full face scan per step;
   the `firstLayer` copy and the count are bitSet operations, and
   `boundFlux` iterates the set bits), and the first call, which plicRDF
   shares, loses the same scan.

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
