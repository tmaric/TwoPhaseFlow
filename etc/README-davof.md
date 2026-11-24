# The `davof-isoAdvector` branch

Reference build of TwoPhaseFlow for the DAVOF (dual area/volume-of-fluid) improvement of
isoAdvector and its verification; the benchmarks live in the sibling `leia` clone
(`<workspace>/leia`). The package is defined in the next section; the branch `davof-swept`
holds the full method development (the swept-face advection, the curved dual advection, the
volume rules) as history.

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
`PATH` and `LD_LIBRARY_PATH`. Nothing is written to, or deleted from, the account default.
`src/VoF/Make/files` lists the davof sources; no third-party library is needed.

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


## The package (branch `davof-isoAdvector`, 2026-11-06)

Tomislav's decision of 2026-11-06: the method (the swept-face advection, the paraboloid
dual advection, the volume rules, the swept geometry) is ported into leia (branch
`davof-port`); this fork keeps ONLY what the isoAdvector improvement needs and runs the
verification that demonstrates it. The branch starts from `davof-swept` at `c5c00c5` and
removes the rest; the diff against upstream `origin/master` (`ba10702`) is the package below
plus the `inline namespace twoPhaseFlow` fix of `feature/parasitic-currents` ("What it is").

(a) The DAVOF reconstruction for isoAdvector (the plane-clipped loop):

- `src/VoF/reconstructionSchemes/plicSchemes/davof/`: the `davof` scheme (the normal from
  the face fractions by the Gauss identity, the plane positioned by alpha; keys
  `wispTol`, `surfCellTol`, `isoFaceTol`, `positioning` (default `ccs`),
  `positioningTol`, `positioningMaxIter`, `normalModel`, `applyTo`, `sliverTol`,
  `maxAngle`, `exactNormal`, `weight`; the per-cell diagnostics `areaNormalDavof`,
  `sliverFlag`, `pPlane`, `planeVolDefect`, `planeAreaDefect`, `positioningResidual`).
- `src/VoF/davof/normalModels/`: the normal models `cellLocal` (the identity normal),
  `youngs`, `plicRDF` (the control), `exact` (the control from the exact state), `ringArea`
  (the stable model: the area-weighted mean over the vertex ring) and `ringFit`
  ("Normal models" below).
- `src/VoF/davof/davofCcsPositioning.H`: the CCS positioning (b).
- `src/VoF/davof/davofInit.{H,C}` with `davofSimplexGeometry.H`,
  `davofQuadraticFaceGeometry.H`, `davofSphereGeometry.H`: the DAVOF state from an implicit
  function (`detrixheAslam`, `quadraticFaces`, `exactSphere`, `exactCylinder`), and `davofOrbit.H` (the
  translated/rotated exact state of the function objects).
- `src/VoF/advectionSchemes/isoAdvection/`: the face update of isoAdvector's own plane
  sweep ("The face update" below): `storeAreaFractions` (the end-of-step submerged face
  areas as `alphaf.<group>`), the second-ring sweep, `planeSpeaks`/`planeFlux`, the
  content bound (`contentBound`), the second bounding layer (`boundSecondLayer`), the
  defaults `defaultAreaFractions`, the processor sync of `alphaf` with `dVf`, the wet-face
  list handed to the scheme; `cutFaceAdvect::subFaceAreaEnd()` and `faceSwept()`.
  Keys of the alpha dictionary: `storeAreaFractions`, `boundSecondLayer`, `contentBound`,
  `debugCell`.

(b) The CCS positioning for every plane-positioning reconstruction scheme
(`src/VoF/surfaceIterators/surfaceIteratorPLIC.{H,C}`, 2026-11-06): `plicRDF` and
`gradAlpha` position their planes with `surfaceIteratorPLIC::vofCutCell`; it now reads
`positioning vofCutCell | ccs;` (default `vofCutCell`, the code unchanged), `positioningTol`
(1e-14 of V_c) and `positioningMaxIter` (30) from the scheme's dictionary
(`readPositioning`). With `ccs` the plane is positioned by `davof::ccsPositionPlane` on the
same `cutCellPLIC` cut (the cubic Hermite interpolant through the bracket ends with the
exact slopes dV/dp = -A, bisection-safeguarded; T. Maric, J. Comput. Phys.: X 11 (2021)
100093, doi:10.1016/j.jcpx.2021.100093, arXiv:2012.02227), started from the linear estimate
between the vertex extrema; the full/empty returns of `vofCutCell` (alpha within `isoFaceTol`
of 0 or 1) are the same on both paths. One Info line per reconstruction with `ccs`
(`surfaceIteratorPLIC: positioning ccs: N positionings, cuts per cell ...`); none with
`vofCutCell`, so the default logs are unchanged. `isoAlpha` positions an iso-value of the
point-interpolated alpha (`surfaceIteratorIso`), not a plane along a normal: the CCS step
needs the slope of the cut volume, which `cutCellIso` does not provide, so `isoAlpha` is
unchanged (the third reference row of the verification). `davof` keeps its own
`positioning` key (default `ccs`).

(c) The function objects `areaFractionError`, `volumeFractionError` (`exactState
exactSphere`), `reconstructionError` (`src/postProcessing/`), the apps `advectorVoF`
(`setFlow` at the half step, the abort guard), `initDavofFields`, `planeVolumeConsistency`,
`staticFaceFractions`, `initShearFlow3D` (the 3D shear flow of Liovic et al. 2006 at t = 0,
U and phi, for the verification's `isoShear3D` case of leia; the time factor is
advectorVoF's `spirallingFlow`), `initVortex2D` (the 2D single vortex of Rider & Kothe
1998 at t = 0 in the x-z plane of a one-cell-thick mesh, phi the exact line integral of
the stream function; the `isoShear2D` case of leia, 2026-11-07), the unit test
`testCutFaceAdvectEndArea`, and the
benchmark dictionaries `run/benchmark/{advection,reconstruction}`; the liquid jet in
cross flow `run/rhoVoF/jetCrossFlow` (its README).

Removed from this branch (history on `davof-swept`): `isoAdvectionSwept.C` and every
`swept*` key, `isoAdvectionCurved.C` (the paraboloid dual advection, the keys
`clippingSurface`, `fluxIntegration`, `fit*`, `anchor*`, `curved*`, `tensorCentroid`,
`planeSource`, `defectLevels`, `writeCurvedFields`), the volume rules
(`src/VoF/davof/volumeRules/`, key `volumeRule`, the vendored overlap library and the Eigen
dependency), the swept geometry (`davofSweptPieces.H`, `davofSweptRemote.H`,
`davofVelocityAt.H`, `davofCurvedGeometry.H`), the face-clipping models
(`src/VoF/davof/faceClipping/`, key `faceClipping`; `planeSagitta` with the curvature
tensor of the `davof` scheme, keys `tensor*`, `writeCurvatureTensor`), the vector face
state of `davof-warped` (`faceState vector`, `wettedAreaDavof`) and `apps/test/testSweptPieces`.
A dictionary that still writes one of these keys is not rejected (unread keys), except
`clippingSurface`/`faceClipping`, which no longer exist as keys at all; `fluxIntegration`
and `volumeRule` are not read.

Identity (keys off = the same code path): the rows davof-ringArea, davof-exact, plicRDF and
isoAlpha of the preprint's smoke cases reproduce the `davof-swept` binary's results bit for
bit; the check and its result are recorded in leia `docs/davof/REPORT-swept.md` ("The
isoAdvector improvement: the verification in the fork (2026-11-06)").

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
4. Boundary faces of a surface cell: the owner's plane sets `alphaf` on an
   outflow face, whose flux is the plane's as before, and on an inflow
   face (phi_f < 0) of a non-coupled patch, whose flux is the boundary
   condition's: the plane's liquid area at the end of the step on both,
   the plane's inflow being the sweep's. The alpha boundary condition's
   value enters the defaults alone (rule 1, `inflowBoundaryFraction`,
   2026-11-06/07): an inflow face no plane sets whose condition says the
   face is wholly liquid (alpha_b = 1: a liquid inlet, the nozzle of the
   jet in cross flow) holds 1 instead of the owner's alpha; every other
   inflow face keeps the owner's alpha, a condition's value below 1 being
   the entering fluid, not the face's state. On the plane gate's inlet
   (fixedValue 0; the plane crosses the inlet and full cells sit at it)
   the versions that took the condition's value in the sweep moved the
   oblique arm in alpha within 20 steps: by 1.8e-2 when every swept inflow
   face took it (690a023), by 5.7e-6 when the liquid inlet's faces took it
   in the sweep, the second ring and the defaults (4da5799, b09b800: the
   plane's cut of an inlet face the condition calls liquid is below 1);
   the defaults alone leave the eight cases unchanged bit for bit (the
   identity check iso/ident5, job 56013017, 2026-11-07). On a processor
   patch the owner's plane
   sweeps the face and the sync replaces the value.
   Decisions 2 and 3 of Tomislav (2026-11-07; the swept loop's keys of
   davof-swept-enright 1306749, taken over for the plane-clipped loop with the
   same names and semantics, the defaults the former behaviour):
   `sweptBoundedFaceFraction none | scale | recompute` -- none: the bounding
   changes a face's dVf and leaves its face fraction; scale: wherever a
   bounding pass changes a face's dVf by s = dVf_bounded / dVf the stored
   face fraction is scaled by s (clipped to [0, 1]; s = 0 dries the face);
   recompute: the face fraction is the end area on the face of the donor's
   plane (the upwind surface cell's, else the second ring's) moved along its
   normal until its swept liquid volume through the face equals |dVf_bounded|
   (a bisection on cutFaceAdvect::timeIntegratedFaceFlux to 1e-13 h), scale
   where the face has no plane or the volume is not bracketed; the counters
   per step ("bounded face fractions (...): faces, s = 0 on, partial, max
   |1 - s|, recomputed (no plane, not bracketed), positioning evaluations").
   `interfaceCellRule normalOnly | normalAndAlpha | normalAndAlphaReset` (the
   davof reconstruction) -- normalOnly: an interface cell is a cell with
   |m_c| > wispTol V^(2/3) regardless of alpha; normalAndAlpha: also
   surfCellTol < alpha < 1 - surfCellTol, a dry or full cell never gets a
   plane and its face fractions are left; normalAndAlphaReset: the same, and
   the faces of such a cell are set to its state (0 / 1) where they differ,
   the shared faces changing the neighbours' m_c too, no mass moving,
   synchronised across processor faces; per step the cells with |m_c| above
   the threshold and alpha outside (tol, 1 - tol), how many got a plane, the
   reset cells and faces and the largest |dalphaf|. The seed log of the shear
   (planeFragmentLog, 2026-11-07) showed the two mechanisms these address: a
   neighbour's extended plane leaving a fraction on a face of a cell without
   liquid (the rule's case) and the bounding correcting a cell to 0 or 1 while
   its faces keep the plane's values (the bounded face fraction's case).
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

(c) The standing rule on plane positioning (Tomislav, 2026-11-08): "We should never count
on the closed form to position PLIC, only always as input for CCS." The explicit
Gauss-theorem position pPlane_c = (sum_f alpha_f (x_f - x_c).S_f^out - 3 alpha_c V_c)/|m_c|
is exact on a planar state only; on every other state it is the START of the CCS
positioning, never the plane. Audit of this branch at bf70393 (2026-11-08): no path takes
pPlane as the final plane. `davof::positionPlane` positions every plane by alpha, with
`positioning ccs` (the default) by `davof::ccsPositionPlane` started from pPlane_c, with
`positioning vofCutCell` by the bisection of `surfaceIteratorPLIC::vofCutCell` (no closed
form); the normal models receive the face-fraction normal m_c (`select(t, m, ...)`), and
the models that need positioned planes (`exact`, `ringFit`) get planes positioned by alpha
in the first pass and re-positioned along the selected normal; the explicit plane is cut
once more for the realizability diagnostics only (`planeVolDefect`, `planeAreaDefect`,
the field `pPlane`); the second-ring sweep of isoAdvection and the `interface` function
object take `surf.normal()`/`surf.centre()`, the positioned plane; `plicRDF` and
`gradAlpha` position by `vofCutCell` or `ccs`; `apps/benchmark/planeVolumeConsistency`
evaluates both positionings by design; `initDavofFields` uses the closed forms of the exact
shapes for the INITIAL alpha and alphaf only (not a PLIC positioning). Any new path that
forms a plane follows the same rule: the closed form may start CCS, never end it.

## The defaults are the PURE set (2026-11-17)

Since 2026-11-17 the code defaults of the davof keys are the PURE set (Tomislav's decision of 2026-11-17): `davofFluxExtension none; davofFaceFractionSource sweepOnly; boundSecondLayer false; sweptBoundedFaceFraction scale; interfaceCellRule alphaAndRingSum; normalModel ringSum; wispTol 0; davofProcFaceRule symmetric;`, the key set of the PURE rows of 2026-11-09 (the `P20o_PURE_*` and `R32_PURE_*` cases of the identity runs). Every key stays selectable. The default gives up the reproduction of the PLIC plane to round-off: with PURE the face fraction on the plane differs from the plane's own cut by 1e-4 to 1e-2 (the range measured on the plane gate of 2026-11-09), and a plane advected by the flux is not reproduced to round-off. The FORMER set is the plane-exact reference and is reached by writing its keys explicitly: `davofFluxExtension secondRing; davofFaceFractionSource secondRing; boundSecondLayer true; sweptBoundedFaceFraction none; interfaceCellRule normalOnly; normalModel cellLocal; wispTol 1e-3; davofProcFaceRule former;`. The old default of `normalModel` was the empty word, which selects `cellLocal` unless the deprecated keys `sliverNormal` or `normalSource` map it onto another model; `normalModel cellLocal` is the same model without that mapping, and the mapping is still reached by `normalModel "";`. The model `ringSum` requires `interfaceCellRule alphaAndRingSum`, so the two keys change together. Only the code defaults change: no algorithm differs, and the post-processing `areaFractionError` keeps its own `wispTol` default of 1e-3.

### Pre-registration of the checks (2026-11-17, before the runs)

Build: this branch `davof-isoAdvector-defaults` (off `davof-isoAdvector` 9abe8c0), new build in `platforms/` of the clone `agents/forkdefault/TwoPhaseFlow`. Old build: `agents/swept/iso/TwoPhaseFlow-h` (9abe8c0, the head of `davof-isoAdvector`, former defaults). Cases: rotation N 32 Co 0.5 (R32, 20 steps) on 1 and 8 ranks and the oblique plane N 20 Co 0.5 (P20o, 20 steps) on 1 rank, as in `agents/swept/iso/ident-entering-sym`, whose `fvSolution` is the template (only the davof key line differs between runs). Compared: every line of the solver log without the clock and execution-time lines, to the digit (`diff`).

- C1: the new build with NO davof key in `fvSolution` gives the same log lines as the same case with the PURE keys written explicitly, on the old build (the stored `*_PURE_*` runs of `ident-entering-sym`, same build and cases) and on the new build (run for R32 on 8 ranks), on all three cases.
- C2: the new build with the explicit FORMER key set gives the same log lines as the old build with no davof key (its former defaults), on all three cases.
- C3: on the 1-rank cases the new build with PURE and `davofProcFaceRule former` gives the same log lines as PURE with `symmetric` (the rule acts on processor faces only); if not, report the difference and stop.
- Runs (12): new build no keys x3, FORMER explicit x3, PURE with rule former x2 (1 rank), PURE explicit x1 (R32, 8 ranks); old build no keys x3. The old build's PURE reference is the stored run of 2026-11-09 (same executable and keys).

### Results of the checks (2026-11-17; build 56534172, runs 56575421, `agents/forkdefault/runs`)

All 12 runs ended with "End" after 20 steps. Compared: the solver log without the clock, execution-time and header lines, `diff`.
- C1 holds: the new build without davof keys gives the same log lines as the old build (9abe8c0) with the PURE keys written explicitly (the stored runs of 2026-11-09) on R32 on 1 rank, P20o on 1 rank and R32 on 8 ranks (0 differing lines on 1 rank; on 8 ranks one header line, the MPI session name `(mpsd0020 8)` against `(mpsd0001 8)`), and the same lines as the new build with the PURE keys explicit on R32 on 8 ranks (0 differing lines). The explicit FORMER keys change 516 of the lines of R32 on 1 rank, so the comparison discriminates.
- C2 holds: the new build with the explicit FORMER key set gives 0 differing lines against the old build with no davof keys, on all three cases.
- C3 holds: on the 1-rank cases R32 and P20o, PURE with `davofProcFaceRule former` and with `symmetric` differ in one line only, the log line that prints the key's own value.
- The build job 56534172 returned rc 1 only from the check `libVoF lacks sweptFaceCells` of the branch `davof`'s build script, a string that does not exist on this branch; the compile and link had no error and every binary exists.

## The ring-sum inconsistency of the ENTERING variant: location (2026-11-23)

Task of Tomislav (2026-11-23, "Reopen the fork track, start with item 2"); branch `davof-isoAdvector-ringsum` off `davof-isoAdvector` ec4e064; workspace
`agents/forkdefault`; tables `etc/ringsum_tables.md`. The isoAdvector-based rows are run at OpenFOAM Co 0.2 and 0.5 only (never 0.95 or 1; standing rule of 2026-11-23).

The record (article FINDINGS 2026-11-09; fork b839296, e854f6e, 9abe8c0): plane translation, normal arm (U parallel to n = (0.3 0.5 0.8124)), Co 0.5, the ENTERING key
set (`davofFaceFractionSource enteringFaces` with its guard): step 1 exact except 3 sliver faces at 1e-5; at step 2 faces off by 6e-3 are written by the interface
cells' own downwind sweeps although the planes are accurate (normal 8.6e-7, centre 1.6e-6), and the reconstruction's own check gives max |A n - m|/V^(2/3) = 0.029
(FORMER 5e-14, PURE 0.16) while the written face fractions are exact; interface-cell counts 176 (ENTERING) against 173 (FORMER); read then: "the ring sum reads data
that differ from the written fields"; never located.

Step 0 (reproduction, the keys written explicitly, job 57087084; the case of the stored 2026-11-09 run `swept/iso/gate-steps/N20_ENTERING_normal_co0.5` with
n_cells, deltaT and endTime changed): the check max |A n - m|/V^(2/3) at the reconstruction after step 1 is 0.02871 at N 20 (stored run: 0.0287146920884307, the
same digits), 0.02577 at N 40, 0.02696 at N 80; 176 interface cells at N 20 after step 1 as stored.

### Pre-registration (before step 1 runs; the three alternatives as given)

- PI1: the 0.029 sits in at most 10 cells, and in each of them m_k "as the code computes it" differs from m_k "from the in-memory alpha_f over all faces" by the mismatch
  itself: a stale or incomplete wet-face list (or cached copy), built before the last modification of alpha_f in the step. Falsifier: the two agree to 1e-12 in those cells.
- PI2 (if PI1 is false): the in-memory field at the check differs from the written field: an order-of-operations fault. Falsifier: they agree to 1e-12.
- PI3 (if both are false): the mismatch is in cells that the active rule admits and FORMER's rule excludes (the 176 against 173), i.e. the check measures sliver cells
  with a poor ring normal, and no datum is read wrongly.
- Report which of the three the data show; the verdict is Tomislav's.

One pre-authorised repair, and only this one: IF PI1 holds in the literal form "a list or cache of alpha_f is not rebuilt after the last modification of alpha_f in the
step", rebuild it at that point (one commit, no key, the line named). Then (i) step 0 again, the check and L1(alpha_f); (ii) the plane ladder N 20/40/80/160, both arms,
Co 0.2 and 0.5, rows FORMER, PURE, ONERING, ENTERING, all keys explicit, before and after the repair in one table (L1(alpha_f) max over steps, the check's max, min alpha,
max alpha - 1); FORMER predicted identical to the digit, PURE to change by less than 10 %, ENTERING normal Co 0.5 to fall below 1e-6. If PURE changes at all, say so
prominently (every stored PURE row would need recomputation: Tomislav's decision). Any other cause than this literal one: no repair, stop and report with the dump of the
worst cell.

### The diagnostics (step 1; print only)

Key `ringSumDiag N` of the davof dictionary (default 0 = off): at every reconstruction of the time indices 1..N, where the check is evaluated (davof.C, after the
normals are selected), the cells with the ten largest |A n - m|/V^(2/3) are printed as `RSDIAG` lines: (a) cell, centre, alpha, interface status under the active rule and under
normalOnly/wispTol 1e-3, |m_k|/V^(2/3); (b) m_k as the code computes it (areaNormal_, through the wet-face list), from the in-memory alpha_f over all faces of the cell (no list),
from the alpha_f field as written to the state's time directory, and from the exact face fractions, and the ring sum with the ring cells; (c) per face: memory, written and
exact alpha_f, membership of the wet-face list, the alphafSource path id (needs `davofFaceSourceLog true`), the sign of phi; (d) n used, n exact, A of the explicit plane's cut,
the CCS residual. Global counters on the `RSDIAG time index` line. The exact state comes from the sub-dictionary `ringSumDiagPlane { origin; normal; U; liquidSide; }`.

### Results (2026-11-23; tables in `etc/ringsum_tables.md`; jobs 57087084 step 0, 57088069/57088093/57088234/57088518/57088519 diagnostics, 57088709 to 57088714 and 57088797/57088798 the key variants)

- Step 0 reproduces the 2026-11-09 record at N 20 to every digit (check 0.0287146920884307, 176 interface cells, L1/Linf(alpha_f) per step); N 40: 0.02577, N 80: 0.02696.
- The check is evaluated in `davof::reconstruct` (davof.C, in the loop over the interface cells: the explicit plane p = (face moment - 3 alpha_k V)/|m_k| along n = m_k/|m_k|, its cut area against m_k; ec4e064 lines 857-874). It does not involve the ring sum: m_k is the cell's own sum, the ring-sum normal is selected after it.
- PI1 is false: in each of the ten worst cells m_k as the code computes it (through the wet-face list) and m_k from the in-memory alpha_f over all faces of the cell agree to 1e-16 (|diff|/V^(2/3) between 1.6e-17 and 1.2e-16); globally, faces with alpha_f != 0 not in the list: 0, faces in the list with alpha_f == 0: 0, at every reconstruction of the time indices 1 to 3. The 0.029 is not in at most 10 cells: 27 of the 176 cells have a check above 1e-3 after step 1 (36 above 1e-6) at N 20; 49 of 143 after step 2.
- PI2 is false: max |alpha_f in memory - alpha_f written to the time directory| = 5.6e-17 at every reconstruction; m_k from the written field equals m_k from memory (|diff| 0 in the printed cells).
- PI3 is false in its literal form: the ten worst cells are interface cells under both rules (alpha_k 0.84 to 0.994); the three cells that only the active rule admits (176 against 173) have alpha_k of 1.6e-6 and 0.999996 and are not among them.
- What the data show instead: the face fractions are plane-exact (max |alpha_f - exact| 1.0e-5 after step 1, from the 3 sliver faces), m_k equals the area vector of the exact plane's cut of the cell (|A_exact| = |m_k| to 1e-12), the planes are accurate (|1 - n.n_exact| below 2e-13 after step 1), but the cell's volume fraction alpha_k is not the exact plane's: alpha_k - alpha(exact plane) = +1.0e-3 (cell 6758: 0.98978 against 0.98877), -1.0e-3, 8.4e-4, ... in the cells of the largest check (33 cells above 1e-6, 2 above 1e-3, at N 20 after step 1; 0.0354 at most after step 2). The explicit plane takes its position from alpha_k, so with exact faces and a wrong alpha_k it cuts an area 17 % off |m_k| in the corner cells: the check. After step 2 the displaced planes (the centre error of 1.6e-6) set the faces by their own downwind sweeps (source id 1) off by 5e-3 to 6e-3: the 6e-3 of the record.
- Where alpha_k departs: single-key variants at N 20 (table): FORMER's keys give alpha_k - alpha_plane below 7e-13 and the check 9e-13 to 8e-12; PURE gives alpha_k off by 1; the ENTERING set keeps its 1.0e-3 when davofFaceFractionSource, boundSecondLayer, sweptBoundedFaceFraction or davofFluxExtension is changed alone, and when interfaceCellRule/normalModel/wispTol go to FORMER's values (check 0.032); with davofFluxExtension secondRing AND davofFaceFractionSource secondRing together (the second-ring sweep sets the flux of the faces beyond the donor), with the ring-sum normal and the alphaAndRingSum rule unchanged, alpha_k is plane-exact (2.4e-13) and the check is 2.5e-12. So the ring sum and the normal model are not the cause; the entering sweep sets alpha_f only ("the flux stays": the volume flux of those faces is the upwind default) while the face fractions come from the plane, and the cell volume fraction advanced with that flux is off by 1e-3 where the plane enters a cell through a face whose upwind cell is not a surface cell.
- No repair was made: the literal condition of the pre-authorised repair (a list or cache of alpha_f not rebuilt after the last modification) does not hold; no ladder was run.

Order within one time step (ec4e064; `isoAdvectionTemplates.C` `advect` unless marked): 484 `surf_->reconstruct()` (davof.C `formAreaNormal` through the wet list from the previous step's `publishWetFaces`, the ring sum, the explicit plane and the check, the normals and the CCS positioning); 503 `timeIntegratedFlux()` (isoAdvection.C 333-745: the sweeps writing alpha_f through `storeAreaFraction`: downwind faces of a surface cell id 1 at 517, upwind faces id 2 at 547, boundary faces ids 3 and 4 at 638 and 670, `sweepSecondRing` at 687 (key source secondRing), `boundCellContent` at 694, `syncProcPatches` at 699); 542 the alpha update; 557 `limitFluxes` (the bounding and, with sweptBoundedFaceFraction scale, the scaling of alpha_f, templates 242 and 277, isoAdvection.C `boundedFaceFraction` 2019); 577 `applyBruteForceBounding`; 600 `sweepEnteringCells` (isoAdvection.C 1035, with `syncEnteringFaces` 1216 for the processor faces); 603 `defaultAreaFractions` (isoAdvection.C 863); 610 `guardFaceFractions` (isoAdvection.C 1285); 623 `publishWetFaces` (isoAdvection.C 983; the bitset `wetFace_` is kept at every writing site: isoAdvection.C 796, 831, 1274, 1322, 1445, 1938, 2117); the solver then writes the fields. The reconstruction of the next step (484) reads the published list.

## The trace of the flux: alpha_k against the exact moving plane (2026-11-23, second task; print only, no repair)

Task of Tomislav (2026-11-23), same branch `davof-isoAdvector-ringsum`; same rules (Co 0.2 and 0.5 only; no repair of any kind in this step).

HYPOTHESIS (Tomislav's, verbatim): the departure of alpha_k from the exact plane is a property of isoAdvector's FLUX, not of DAVOF: isoAdvector integrates the wetted area only on the downwind faces of cells that hold a plane at t^n; a face whose upwind cell holds no plane at t^n (full or empty) takes alpha_up phi_f dt even when the plane, coming from a neighbouring cell, sweeps over part of that face during the step. The flux error of such a face is the corner volume the plane cuts from its swept prism.

- PT1: every face with |dVf - dVf_exact| > 1e-12 V_k in step 1 (ENTERING and PURE alike) has an upwind cell WITHOUT a plane at t^n; no downwind face of an interface cell has such an error. Falsifier: one downwind face of an interface cell above 1e-10 V_k.
- PT2: summed over the faces of a cell, these face errors reproduce the cell's alpha_k deviation from the exact plane to 1e-12 (so nothing else contributes: not the bounding, not the limiter).
- PT3: plicRDF + isoAdvector and isoAlpha + isoAdvector on the SAME cases show the same faces with the same flux errors in step 1 (identical initial alpha; their planes are exact on a plane to their tolerance), and alpha_k deviations of the same size over the run as PURE (within a factor 2 in the max norm and in L1).
- PT4: at Co 0.2 the number of such faces and the max alpha_k deviation are smaller than at Co 0.5 but not zero on the normal arm; on the oblique arm (U = (1 0 0)) report what is found (no prediction).

MEASURE. A: step 1, N 20 and 40, normal arm Co 0.5, rows ENTERING and PURE: per face with error above 1e-12 V_k: label, dVf, exact, error/V_k, upwind cell (label, alpha^n, plane at t^n yes/no), downwind cell, the alphafSource id; counts by class; per cell the sum of its faces' errors against its alpha_k deviation (PT2), top 10. B: the per-step counters (faces above 1e-12, max and L1 of the face flux error, max and L1 of the alpha_k deviation, the explicit-plane check where it exists) over the whole run for the rows FORMER, PURE, ONERING, ENTERING, plicRDF, isoAlpha on N 20/40/80, both arms, Co 0.2 and 0.5; for plicRDF and isoAlpha also the normal error at the end.

HOW the exact flux is computed (key `ringSumFluxTrace N`, default 0 = off; counters at every step, details for the time indices 1..N; `isoAdvection::fluxTraceReport`, isoAdvection.C, called at the end of the flux step after the bounding and the clipping; the exact state from `ringSumDiagPlane {origin; normal; U; liquidSide}`): the exact flux of face f over [t^n, t^{n+1}] is phi_f times the time integral of the exact wetted fraction of the face for the plane x0(t) = origin + U t. The fraction of a planar face is piecewise quadratic in t between the times at which the plane passes a vertex of the face, so Simpson's rule on the intervals between those times is exact (a face on one side of the plane at both ends of the step has fraction 0 or 1 throughout). dVf is compared before the limiter (`dVfRaw`, copied after `timeIntegratedFlux`) and after it. The exact volume fraction of a cell is the exact cut of the cell by the plane (cutCellPLIC, orientation fixed against alpha^0 at the first call; cells with all vertices on one side take 0 or 1). A plane exists in cell c at t^n if `surf_->interfaceCell()[c]` and its normal is not zero (the test timeIntegratedFlux uses). The cell sum E = - sum_f sigma_f (dVf_f - exact_f)/V_c (sigma = +1 for the owner) is the alpha_k deviation the face errors give when alpha^n is exact; it is printed next to the measured alpha_k - alpha_plane(t^{n+1}) (after the limiter and the clipping) so that the limiter and the clipping show as the difference.

### Results of the flux trace (2026-11-23; tables `etc/ringsum_tables.md`, jobs 57089559 (test) and 57089578 to 57089606 (the 24 ladder jobs); build 57089456; no repair)

Rows x arms x Co x N: 72 runs (6 rows, normal and oblique arm, OF Co 0.2 and 0.5, N 20/40/80); every run ended with "End". N 160 was not run. The initial state of every row is the initDavofFields state (the comparators too), not the setAlphaField state of the earlier studies.

- Step 1, normal arm, Co 0.5, ENTERING, PURE and ONERING (the same numbers in the three rows): N 20: 21 faces with |dVf - exact| > 1e-12 V_k (19 internal, 2 boundary); N 40: 89 (84 + 5); N 80: 386 (376 + 10). In every one of them the upwind cell has no plane at t^n (alpha^n = 1 to 1e-8) and the downwind cell has none; the largest error is 1.01e-3 V_k (N 20), 8.6e-4 (N 40), 9.2e-4 (N 80); the faces whose upwind cell holds a plane have errors of at most 4.8e-15 V_k (N 20) and 1.5e-13 V_k (N 40). FORMER: no face above 1e-12 V_k in step 1 (max 1.7e-14).
  Example, N 20, face 19614 (upwind cell 6757, alpha^n = 0.99999999236, no plane; downwind cell 6758): dVf = 1.2000e-12 m^3 = alpha_up phi_f dt, exact 1.1919e-12 m^3, error/V_k = 1.011e-3.
- PT2 at step 1: in the ten cells with the largest alpha_k deviation D (cells of alpha^n = 1, no plane), the sum of the face errors E = - sum_f sigma (dVf - exact)/V_k equals D to 1.8e-15 (N 20), 7.1e-15 (N 40), 1.5e-14 (N 80); the faces above 1e-12 V_k are the same before and after the limiter (21 and 21 at N 20); alpha^n - alpha_plane(t^n) is below 5e-15; sum |D| V_k = 9.60e-14 m^3 at step 1, N 20, the value of volumeFractionError's E_shape.
- From step 2 on, errors also sit on faces whose upwind cell holds a plane (N 20, ENTERING, step 2: 61 faces up, 54 both, 18 none; step 3: 71, 247, 19): the planes of the interface cells are positioned from the alpha_k that already deviates.
- Comparators (same initial alpha): plicRDF + isoAdvector and isoAlpha + isoAdvector have 530 and 543 faces above 1e-12 V_k in step 1 at N 20, Co 0.2 (normal arm), and 623 and 629 at Co 0.5, of which 384 have a plane both upstream and downstream; the alpha_k deviation after step 1 is 2.8e-3 (Co 0.2) and 3.4e-3 (Co 0.5) for plicRDF and 4.3e-2 and 9.2e-2 for isoAlpha, against 1.0e-3 for PURE at Co 0.5; alpha_k deviation over the run (normal arm, Co 0.5, N 20/40/80): PURE 6.5e-3/5.8e-3/6.3e-3, plicRDF 3.8e-2/8.6e-2/1.1e-1, isoAlpha 0.19/0.19/0.19 (Co 0.2: PURE 1.2e-4/1.2e-4/1.3e-4, plicRDF 2.8e-3/2.8e-3/2.5e-3, isoAlpha 0.23); normal error of the comparators at the end (LNormalInf_ / LNormal1_): plicRDF 3.5e-10/7.6e-11 (N 20), 5.1e-14 (N 80) at Co 0.2 and 0.031/0.006 (N 20), 0.16/0.013 (N 40), 0.095/0.007 (N 80) at Co 0.5; isoAlpha 0.19/0.076 (N 20, Co 0.2).
- Co 0.2, normal arm: PURE has 3, 0 and 0 faces above 1e-12 V_k in step 1 (N 20/40/80; the 3 at 2.5e-9 V_k), ENTERING the same; the alpha_k deviation over the run is PURE 1.2e-4/1.2e-4/1.3e-4, ONERING the same, ENTERING 1.1e-7/2.7e-14/6.5e-14, FORMER 2.0e-14/2.1e-13/3.6e-10. Oblique arm, Co 0.2: step 1 has 1, 0 and 0 faces (PURE, ONERING, ENTERING); deviation over the run PURE 4.7e-2/5.8e-2/7.9e-2, ONERING 3.9e-6/5.4e-6/5.7e-6, ENTERING 4.7e-8/7.8e-8/7.8e-8, FORMER 4.2e-14/3.3e-12/1.4e-9; Co 0.5: PURE 2.8e-2/4.3e-2/5.7e-2, ONERING 9.5e-5/1.2e-4/1.7e-4, ENTERING 6.8e-8/8.0e-14/3.2e-13, FORMER 4.4e-14/3.3e-12/5.0e-9.
- Unexpected: FORMER on the normal arm at Co 0.5 is exact in step 1 (1.9e-14) and then grows by a factor of 20 to 40 per step: N 40, alpha_k deviation 5e-13 (step 2), 1.8e-11, 5.5e-10, 3.0e-8, 9.3e-7, 3.7e-5, 1.4e-3 (step 8), 0.86 at the end (16 steps); N 80: 1 at the end (32 steps), N 20: 8.6e-9 after 8 steps; the explicit-plane check grows with it (1.5e-13 at the first reconstruction, 4e-8 after step 6 at N 40). At Co 0.2 the growth is slower (N 80, 80 steps: 3.6e-10).
- The numbers of the plane-ladder rows and the comparators by row, arm, Co and N are in the table B of `etc/ringsum_tables.md`.

## The Courant number of the plane ladder; the normal arm at OpenFOAM Co 0.2 and 0.5 (2026-11-23, third task; print only, no repair, no new key)

Task of Tomislav (2026-11-23). Step 1, read-only: the labelled Co of the ladder is |U| dt/h (the step counts are the same on both arms); OpenFOAM's Courant number is (|u_x| + |u_y| + |u_z|) dt/h. Measured from the solver's own `Courant Number mean/max` lines, N 20 (the stored 2026-11-09 case has the same deltaT = 1e-4 at N 20):

| arm | label | deltaT (N 20) | solver Courant Number mean / max | ratio to the label |
|---|---|---|---|---|
| normal, U = (0.3 0.5 0.8124) | 0.5 | 1e-4 | 0.8062 / 0.8062 | 1.6124 |
| normal | 0.2 | 4e-5 | 0.32248 / 0.32248 | 1.6124 |
| oblique, U = (1 0 0) | 0.5 | 1e-4 | 0.5 / 0.5 | 1.0 |
| oblique | 0.2 | 4e-5 | 0.2 / 0.2 | 1.0 |

So the normal-arm rows labelled "Co 0.5" and "Co 0.2" of the ladder (and the 2026-11-09 record) are OpenFOAM Co 0.806 and 0.322. Step 2 reruns the normal arm at OpenFOAM Co 0.2 and 0.5, i.e. |U| dt/h = 0.2/1.6124 = 0.124039 and 0.5/1.6124 = 0.310097.

Pre-registration of step 2 (verbatim from the task, before the runs):
- PE1: at OpenFOAM Co 0.2 and 0.5 on the normal arm ENTERING has 0 step-1 faces above 1e-12 V_k at N >= 40 and a max alpha_k deviation over the run <= 1e-10 at N >= 40 (N 20: <= 1e-6, the guard-dried slivers). Falsifier: one N >= 40 case above 1e-8.
- PE1b: FORMER stays below 1e-8 at these Courant numbers at every N (its growth at 0.806 is outside the range).
- PE1c: PURE's L1 deviation is below plicRDF's on the normal arm and above it on the oblique arm (as at the old labels).

Definitions of the normal norms in the tables (reconstructionError.C): LNormalDiff1 and LNormalDiffInf are the L1 (mean over the interface cells) and Linf of |n - n_e| at the last step; LNormal1_ is the time average over the run of the per-step mean of 1 - n.n_e and LNormalInf_ the maximum over the run of its maximum (so the 'LNormalInf_ / LNormal1_' of the earlier tables are run quantities, not end-of-run values).

### Results of step 2 (2026-11-23; tables in `etc/ringsum_tables.md`, section "The plane ladder with the Courant numbers labelled"; jobs 57090102 to 57090113; 36 runs, all ended with "End")

- The normal-arm rows of the earlier table B are relabelled: |U| dt/h = 0.2 and 0.5 are OpenFOAM Co 0.322 and 0.806 (0.806 outside the range); the rows rerun at OpenFOAM Co 0.2 and 0.5 (|U| dt/h = 0.124039 and 0.310097; the solver prints Courant Number 0.19999999988 at N 40, OpenFOAM Co 0.2) are in the new table next to them. The oblique arm's labels are OpenFOAM Co.
- PE1 is false. ENTERING, normal arm, max alpha_k deviation over the run at N 20/40/80: OpenFOAM Co 0.2: 7.1e-8, 5.0e-8, 7.8e-8; OpenFOAM Co 0.5: 2.6e-7, 2.7e-7, 2.8e-7 (the falsifier, one N >= 40 case above 1e-8, is met by all four N >= 40 cases). Step-1 faces above 1e-12 V_k: Co 0.2: 3, 0, 0; Co 0.5: 3, 4, 20 (N >= 40 at Co 0.5 not zero). At OpenFOAM Co 0.2, N 40, ENTERING is exact (1.2e-14, 1.5e-14) at steps 1 and 2; at step 3 the deviation is 3.9e-8 with 33 faces above 1e-12 V_k, and stays at 2e-8 to 5e-8 (1500 to 1700 such faces per step) through step 14 (the steps inspected; the run has 64 steps, its maximum is 5.0e-8); at the old label (OpenFOAM Co 0.322) the same N 40 case stays at 2e-14 for all 40 steps.
- PE1b: FORMER is below 1e-8 at OpenFOAM Co 0.2 (1.4e-14, 1.4e-10, 1.2e-9 at N 20/40/80) and at OpenFOAM Co 0.5 for N 20 (8.3e-9); at OpenFOAM Co 0.5 it is 0.955 at N 40 and 1.0 at N 80 (the growth seen at 0.806 occurs at 0.5 as well), so PE1b is false at OpenFOAM Co 0.5, N >= 40.
- PE1c: on the normal arm PURE's L1 deviation (sum |D| V / V_domain) at the last step is above plicRDF's at OpenFOAM Co 0.2 and 0.5 and at 0.322 (1.4e-8 against 8.1e-9 at Co 0.2 N 20; 5.2e-8 against 2.4e-10 at Co 0.5 N 80), below plicRDF's at 0.806; the maximum over the steps is below plicRDF's at Co 0.2 and 0.322 and in 2 of 3 N at 0.5. On the oblique arm PURE is above plicRDF's in every case (last step and maximum), as predicted. So PE1c holds on the oblique arm and at 0.806 on the normal arm, and is false at OpenFOAM Co 0.2 and 0.5 on the normal arm for the last-step L1.
- Step 3 (the curved standing matrix) was not run: PE1 fails, as agreed.

POST HOC (2026-11-23, main session): "PE1's falsifier (1e-8) was set by me without regard to the scheme's own interface band surfCellTol = 1e-8; ENTERING's measured deviation on the normal arm, 5e-8 to 2.8e-7 at OpenFOAM Co 0.2 and 0.5, falsifies PE1 as written. It is four orders below plicRDF's maximum deviation (2.5e-3 to 2.8e-3) and three below PURE's. The gate for the curved matrix is amended post hoc to: ENTERING's max alpha_k deviation on the plane <= 1e-6 at OpenFOAM Co 0.2 and 0.5, both arms, N 20 to 80, which the data meet. PE1 stays recorded as falsified; PE1b is falsified at OpenFOAM Co 0.5 (FORMER 0.955 and 1.0 at N 40, 80); PE1c is falsified on the normal arm in the last-step L1."

### Step A: location of the step-3 jump (2026-11-23, fourth task; print only, existing keys only)

Case: ENTERING, normal arm, OpenFOAM Co 0.2 (|U| dt/h = 0.124039), N 40: exact (1.2e-14, 1.5e-14) at steps 1 and 2, 3.9e-8 with 33 faces above 1e-12 V_k at step 3.
- HJ1 (recorded before the run): the 33 faces at step 3 have an upwind or downwind cell whose exact alpha lies within the band (0, 1e-8) or (1 - 1e-8, 1) at t^n, i.e. a cell the plane cuts that the rule treats as full or empty (no plane), or a face the guard dried or filled; the deviation per cell is of the order of the band times a geometric factor and does not grow beyond 1e-6 over the run. Falsifier: a face of the 33 with both cells outside the band and planes where the rule says so.
- HJ2: with the existing tolerance key(s) of that band lowered (surfCellTol 1e-8 -> 1e-10 and 1e-12, and the guard's eps if it is an existing key, one at a time and together; existing keys only) the run's maximum deviation falls roughly in proportion. If the scan needs a tolerance that is not an existing key: no key is added, report.
The faces are printed by the flux trace with `ringSumFluxTrace 3` (the line now also gives the exact alpha^n of the upwind and downwind cells and the stored alpha^n of the downwind cell).

### Results of step A: the step-3 jump (ENTERING, normal arm, OpenFOAM Co 0.2, N 40; jobs 57090450 to 57090452; build 57090424)

- HJ1 is falsified as written: of the 33 faces above 1e-12 V_k at step 3 (the largest 3.8e-8 V_k), 30 have both the upwind and the downwind cell holding a plane at t^n, with exact and stored alpha^n in the interior of (0, 1) (0.0089 to 0.9986, the stored value within 1e-14 of the exact one); only 3 faces (136799, 136800 and the boundary face 189958) have cells without a plane, and there the upwind cell has alpha^n = 0.9999999975, inside the band (1 - 1e-8, 1). The error sits on faces whose upwind cell has a plane (the largest, face 127439, 3.77e-8 V_k: upwind cell 43199, alpha^n = 0.1485, plane; downwind cell 44799, plane).
- What the log shows at the step before: at step 2 the guard acts: "seeds before the guard: 13 cells with alpha = 0 and a wet face, 2 cells with alpha = 1 and a face below 1", "guard: 15 faces set to 0 (0 above 1e-12), 6 to 1 (6 above 1e-12)" (at step 1: none). The reconstruction error of the normals is 1.2e-14 after step 1 and 4.6e-8 after step 2, 1.0e-7 after step 3; the explicit-plane defect |alpha_plane - alpha| and the check are 5.5e-14 and 5.4e-13 at the reconstruction after step 1 and 2.5e-6 and 6.2e-6 at the one after step 2. So the guard's face writes of step 2 (the second clause of HJ1, "a face the guard dried/filled") are followed one step later by the flux errors on the faces of the neighbouring interface cells; the mechanism between the two (the ring sum of m over the ring cells that read those faces) was not traced separately.
- HJ2 (existing key surfCellTol, the band of the entering sweep and the guard: the guard's eps is surfCellTol_, no separate key exists, so the three variants "one at a time and together" are one): ENTERING, normal arm, OpenFOAM Co 0.2, N 40, 64 steps:

| surfCellTol | max alpha_k deviation over the run | L1 rel (last / max over steps) | max check | faces above 1e-12 V_k at step 3 | deviation at step 3 |
|---|---|---|---|---|---|
| 1e-8 (standing) | 5.0e-8 | 4.1e-12 / 4.6e-11 | 1.4e-5 | 33 | 3.9e-8 |
| 1e-10 | 1.3e-9 | 5.2e-14 / 4.6e-13 | 8.8e-7 | 0 | 1.4e-10 |
| 1e-12 | 1.4e-10 | 3.6e-15 / 7.1e-15 | 8.8e-7 | 0 | 1.4e-10 |

  The maximum falls by a factor 38 for a factor 100 in the tolerance, and by a further factor 10; the step-3 deviation falls from 3.9e-8 to 1.4e-10 and stays there for 1e-12.

### Pre-registration of step B (verbatim from the task, before the runs): the curved standing matrix

Shear2D N 64/128/256, shear3D N 32/64/128, deformation3D N 32/64/128 (the same ranks as the standing rows: 1/1/8, 1/8/64, 1/8/64), OpenFOAM Co 0.2 and 0.5 only, rows plicRDF, PURE, ONERING, ENTERING, every davof key explicit, davofProcFaceRule symmetric in all davof rows, the fork binary of this branch through the iso-time.sbatch route of the 2026-11-09 matrix (leia ac5fd7f5 / 9527a46f), one repetition per case, a version stamp per case; standing tolerances (surfCellTol 1e-8 etc.).
- PE2: ENTERING completes all 18 cases (falsifier: one abort).
- PE3: E_shape of ENTERING <= PURE in at least 14 of 18 and ENTERING / plicRDF between 0.9 and 1.5.
- PE4: wall clock of ENTERING <= 1.1 x plicRDF at equal ranks in every case.
- PE5: ONERING completes all 18 and lies between PURE and ENTERING in E_shape in at least 12 of 18.
If ENTERING aborts anywhere: no repair; the first abort is reported with the abort line, the step, min/max alpha, and the five cells with the largest violation at the step before; the other cases continue.

### Results of step B: the curved standing matrix (2026-11-23; `etc/entering_matrix_tables.md`; jobs 57090666 to 57090683, 18 jobs of 4 rows, one repetition; fork binary of this branch, build 57090424, stamp per case in `agents/forkdefault/iso/time/*/*_r1/.stamp`)

- PE2 holds: ENTERING completed all 18 cases (no abort). ONERING and PURE and plicRDF also completed all 18.
- PE3 is false in its first part: E_shape of ENTERING is below PURE's in 3 of 18 cases (ENTERING / PURE: 0.9999, 0.9999, 0.9996). In the three-dimensional flows ENTERING and PURE (and ONERING) agree to 1e-3 or better in 10 cases (ratios 0.9996 to 1.0014; shear3D and deformation3D give E_shape 0.00629, 0.00188, 5.0e-4 and 0.0100, 0.0039, 7.6e-4 in all three rows); in shear2D ENTERING is above PURE by 2 % to 25 % (1.019 to 1.251). ENTERING / plicRDF is between 0.9 and 1.5 in 13 of 18 cases: below 0.9 in shear2D N 64 (0.60, 0.73), above 1.5 in shear3D N 128 (1.52, 1.50 at Co 0.2, 0.5; 1.34 to 1.41 at N 32 and 64) and between for the rest (0.88 to 1.10 in shear2D N 128 and 256, 1.10 to 1.32 in deformation3D).
- PE4 holds: wall clock of ENTERING / plicRDF at equal ranks is 0.45 to 0.99 (0.99 in shear3D N 32 at both Co; 0.57 to 0.82 otherwise); one repetition per row, rows run one after the other on one exclusive i01 node.
- PE5 is false in its second part: ONERING completes all 18, but lies between PURE and ENTERING in E_shape in 2 of 18 cases.
- Courant numbers (the solver's own maximum, labels 0.2 / 0.5): shear2D 0.1996 to 0.1998 / 0.4991 to 0.4996; shear3D 0.2165 to 0.2170 / 0.5413 to 0.5426; deformation3D 0.1968 to 0.1983 / 0.4920 to 0.4956. The solver does not print the target in the shear3D rows (+8 % to +9 %) and in deformation3D (-1 % to -1.6 %).
- The tables hold, per row and case: status, E_shape and order, L1(x) and order, |n - n_e| L1 (with order) and Linf, 1 - n.n_e L1 and Linf, mass error, min alpha, max alpha - 1, interface cells at the end, wall clock.

## Follow-ups of 2026-11-24 (Tomislav's order): the 3D shear gap, the sphere along an axis, the lower bound

### Correction of the normal columns (2026-11-24)
The column "abs(n-n_e) L1" of `entering_matrix_tables.md` read 1/N^3 (3.05e-5, 3.81e-6, ...): the column index into `reconstructionError`'s `error.dat` was one too low, so it read `avgVol` (the cell volume) instead of `LNormalDiff1`; the columns "abs(n-n_e) Linf", "1-n.n_e L1" and "1-n.n_e Linf" read `LNormalDiff2`, `LNormalDiff1` squared / 2 and `LNormalDiff2` squared / 2 for the same reason. The columns are recomputed from the right ones (`LNormalDiff1` = L1 of |n - n_e|, `LNormalDiffInf` = Linf, and 1 - n.n_e from `LNormalDiff2` (the RMS) and `LNormalDiffInf` as |n - n_e|^2/2 for unit normals) in the committed matrix table; the same index error was in the comparator columns of the ladder table of the third task (`ringsum_tables.md`, section "The plane ladder with the Courant numbers labelled"), regenerated; the columns "LNormalInf_ / LNormal1_" of the first and second reports were read from the right columns.

### Task 3 (print only): the lower bound of the DAVOF rows
Key `ringSumBoundTrace N` (default 0 = off; `isoAdvection::limitFluxes`, isoAdvectionTemplates.C): per step the min and max - 1 of alpha before the bounding, after each pass and at the end of the step; the cells outside [-1e-12, 1 + 1e-12] after the last pass inside and outside needBounding (surfCells_ plus one layer; the davof rows have boundSecondLayer false); the five worst cells of a step with violation above 1e-6 (alpha before and after, needBounding, plane at t^n, outflow faces with dVf != 0). Case: deformation3D N 64, OpenFOAM Co 0.5 (target label 0.5), rows plicRDF and PURE, serial.
HYPOTHESIS (HB1, verbatim): in PURE the cells left below zero after the last pass lie OUTSIDE needBounding in the majority of steps where min alpha < -1e-6 (the DAVOF interface-cell set, with wispTol 0 and the alphaAndRingSum rule, differs from plicRDF's surfCells_), while in plicRDF they lie inside. Falsifier: the majority inside needBounding in PURE too; then report the pass-by-pass numbers (bounding not converged in three passes, or the violation larger before the bounding).

### Task 1: why the davof rows are above plicRDF on the 3D shear
Controls, shear3D N 32/64/128, OpenFOAM Co 0.2 and 0.5, the ranks of the matrix: (a) plicRDF; (b) plicRDF with positioning ccs; (c) the PURE loop with the normal of plicRDF's normal computation instead of the ring sum, everything else PURE (a control, not a method); (d) PURE. Diagnostics D1: at t = 0 on the exact initial state the normal errors (1 - n.n_e and |n - n_e|, L1 and Linf) of the ringSum normal and of plicRDF's, N 32/64/128 with orders (the error row of a run with deltaT/1000 and one step; the state is the initial one up to 1e-3 of a step); D2: L1 and Linf of alpha_f against the exact sphere's face fractions at t = T for PURE and right after the initialisation, with orders; D3: plicRDF's normal computed on PURE's final alpha (a restart of the plicRDF case from PURE's final fields with deltaT/1000 and one step), its 1 - n.n_e next to plicRDF's own run and PURE's own.
PRE-REGISTERED (verbatim): PN1 statically (D1) the ringSum normal's 1 - n.n_e L1 is within 15 % of plicRDF's at every N. PN2 row (c) has E_shape within 5 % of plicRDF's and of row (b) at every N and Co: the normal source accounts for the whole gap. PN3 (D3) plicRDF's normal on PURE's final alpha has 1 - n.n_e L1 within 20 % of plicRDF's own run: the alpha fields are comparable and the excess normal error of PURE is in the transported face fractions. Falsifiers: PN1 beyond 30 %; PN2 a remaining gap above 15 %; PN3 beyond 50 %. Report each.

### Task 2: does PURE's axis-aligned plane defect show on a curved interface
Sphere translation along a mesh axis: the fork's sphereTranslation3D case (box 2 x 1 x 1, R = 0.15 at (0.5 0.5 0.5), U = (1 0 0), end time 1 = 3.3 diameters, exact initial state by initDavofFields with the exact sphere), N 32/64/128, OpenFOAM Co 0.2 and 0.5 (U = (1 0 0): OpenFOAM Co = |U| dt/h), rows plicRDF, PURE, ONERING, ENTERING.
PRE-REGISTERED (verbatim): PS1 E_shape of PURE is at least 20 % above ENTERING's and ONERING's at N 64 and 128, both Co (the defect shows). PS2 ENTERING's E_shape is <= plicRDF's at every N and Co. PS3 all rows complete with alpha within 1e-6 of [0,1]. Falsifier of PS1: PURE within 10 % of ENTERING everywhere (then the plane defect is immaterial on this curved case).

### Results of the follow-ups of 2026-11-24 (tables: `etc/bounding_tables.md`, `etc/normal_control_tables.md`, `etc/sphere_axis_tables.md`; jobs: bounding 57098598 (a first trace run was cancelled and rerun after the list append of the trace was made linear), D1 57094840/57094849/57094850, controls 57095132 to 57095137, D3 57096995 to 57097000, sphere 57094999 and 57095143 to 57095147; build 57098024)

Task 3 (the lower bound), deformation3D N 64, target Co 0.5, serial, nAlphaBounds 3:
- plicRDF: min alpha over the run is -0.0118 before the bounding and -2.1e-4 after the last pass; PURE: -0.0287 before and -5.3e-3 after (25 times plicRDF's); max alpha - 1: 0.018 / 3.1e-4 (plicRDF), 0.026 / 2.0e-5 (PURE).
- HB1 (the cells below zero after the last pass lie outside needBounding in the majority of the steps with min alpha < -1e-6 in PURE, inside in plicRDF): PURE has 265 such steps; in 200 of them the worst cell is outside needBounding and in 257 more cells above 1e-6 lie outside than inside. plicRDF has 67 such steps; the worst cell is outside in 17 of them and more cells lie outside than inside in 26. The five worst cells of the run in PURE are inside needBounding (3 of 5 are interface cells at t^n, 2 hold no plane at t^n; alpha -5.3e-3, -4.4e-3, -3.9e-3, -3.2e-3, -3.0e-3 after the bounding, with 0 to 3 outflow faces carrying a flux), as are plicRDF's (two cells at 1 + 3.1e-4, one at -2.1e-4, two at -1.6e-4, all interface cells at t^n).
- Pass by pass at the step of PURE's lowest alpha (step 922): before the bounding min -0.0225, after pass 1 -0.00675, after pass 2 -0.00803 (lower than after pass 1), after pass 3 -0.00527; cells outside [-1e-12, 1 + 1e-12] after the last pass: 1007 inside needBounding and 65627 outside. plicRDF's step of lowest alpha (59): -0.00664, -0.00175, -0.000778, -0.000205; 162 inside and 2501 outside. The bounding is not converged after three passes in both rows; in PURE a pass can lower the minimum.

Task 1 (the 3D shear): D1 at t ~ 0 on the exact initial state (one step of deltaT/1000): 1 - n.n_e L1 of the ringSum normal is 7.9e-4, 3.0e-4, 7.4e-5 (N 32/64/128; orders 1.4, 2.0) and plicRDF's 4.2e-4, 1.5e-4, 3.5e-5 (orders 1.5, 2.1): ratios 1.90, 2.06, 2.10; |n - n_e| L1 0.0318, 0.0178, 0.0088 against 0.0207, 0.0123, 0.0060. PN1 is false (beyond 30 %). D2: L1(alpha_f) of PURE against the exact sphere's face fractions 0.0070, 0.0034, 0.0019 right after the initialisation (Linf 0.036, 0.023, 0.017) and 0.26, 0.17, 0.10 (Co 0.2) and 0.24, 0.16, 0.094 (Co 0.5) at t = T (Linf 1), orders 0.6 to 0.8. Row (c) (the PURE loop with plicRDF's normal, the existing model normalModel plicRDF; no model was added): E_shape within 0.1 % of plicRDF's and of row (b) (plicRDF with ccs) at all six cells (ratios 1.000 to 1.001), against PURE's 1.335 to 1.523 times plicRDF's: PN2 holds. D3: plicRDF's normal computed on PURE's final alpha has 1 - n.n_e L1 1.47 to 1.70 times plicRDF's own run (1.63, 1.47, 1.61 at Co 0.2; 1.70, 1.48, 1.48 at Co 0.5), and 1.05 to 1.08 times PURE's own: PN3 is false (beyond 20 %; beyond 50 % in 3 of 6).

Task 2 (the sphere along an axis): all 24 runs completed with alpha within 1e-12 of [0, 1] (PS3 holds); the solver's Courant maximum prints the label. E_shape, N 32/64/128: Co 0.2: plicRDF 3.23e-4, 8.11e-5, 2.16e-5; PURE 2.78e-4, 7.35e-5, 2.12e-5; ONERING 2.78e-4, 7.34e-5, 2.12e-5; ENTERING 2.78e-4, 7.33e-5, 2.12e-5. Co 0.5: plicRDF 4.87e-5, 1.48e-5, 2.69e-6; PURE 6.26e-5, 1.02e-5, 2.73e-6; ONERING 6.01e-5, 9.80e-6, 2.28e-6; ENTERING 5.94e-5, 9.60e-6, 1.93e-6. PS1 (PURE at least 20 % above ENTERING and ONERING at N 64, 128, both Co) is false: PURE / ENTERING is 1.002 and 1.002 at Co 0.2, 1.060 and 1.415 at Co 0.5 (PURE / ONERING 1.001, 1.002, 1.038, 1.198); the falsifier (PURE within 10 % everywhere) is not met either (N 128, Co 0.5). PS2 (ENTERING <= plicRDF at every N and Co) holds in 5 of 6 (ENTERING / plicRDF 0.86, 0.90, 0.98 at Co 0.2; 1.22, 0.65, 0.72 at Co 0.5).

## Branch davof-isoAdvector-normal (off 826d8bb), 2026-11-23: the static normal gate and the bounds control; pre-registration (Tomislav, written before the first run)

Rules: OpenFOAM Co 0.2 and 0.5 only, checked against the solver's own Courant line; every davof key explicit in every case; one binary for all rows of a table; nothing repaired, tuned or added beyond what is specified; every prediction is reported as held or false with the numbers.

### Part N: static normal gate (no transport)

Route: the exact initial state by initDavofFields, one step of deltaT/1000, reconstructionError. States: the 3D sphere (shear3D initial state) at N 32/64/128 at the standard centre and at the centre shifted by (0.013 0.007 0.003); the 2D circle (shear2D initial state) at N 64/128/256/512 at the standard centre and at the centre shifted by (0.013 0.007 0). Rows (the normal only, everything else as PURE): (a) plicRDF; (b) normalModel ringSum; (c) the ring quadric (ringFit on the old davof branch): run if it exists on 826d8bb under any name, otherwise reported and skipped, not ported; (d) the new control model normalModel ringSumRDF with key rdfIterations (1 and 2 are run): step 1 the ring-sum normal in every interface cell, the plane positioned from alpha_k as PURE does; step 2 the reconstructed distance function from these planes on the point-neighbour stencil and the normalised least-squares gradient of it at the interface cells, by the same class and calls plicRDF uses for its update, the plane repositioned from alpha_k; step 2 repeated rdfIterations times; no residual tolerance, no fallback, no other branch. If plicRDF's update cannot be called this way without changing plicRDF's own code path: stop and report what blocks it. Source header cites H. Scheufler, J. Roenby, "Accurate and efficient surface reconstruction from volume fraction data on general meshes", J. Comput. Phys. 383 (2019) 1-23, doi:10.1016/j.jcp.2019.01.009, arXiv:1801.05382 (checked against plicRDF's header). Check first at N 32: rdfIterations 0 reproduces row (b) to the digit. Reported per state, position, N and row: 1 - n.n_e L1 and Linf, |n - n_e| L1 and Linf, orders, the ratio to row (a); the number of RDF iterations plicRDF takes in this step; plicRDF's mean iteration count per step in the logs of the shear3D matrix runs, if they print it.

- PG1: the ring quadric's 1 - n.n_e L1 on the 3D sphere is at most 0.5 times plicRDF's at N 64 and 128, both positions. Falsifier: above 1.0 in any of the four.
- PG2: ringSumRDF with one iteration is within 10 % of plicRDF's 1 - n.n_e L1 at every N, both states, both positions. Falsifier: beyond 25 % anywhere.
- PG3: on the 2D circle the static ring-sum normal is not better than plicRDF's: the ratio of |n - n_e| L1 is at least 1.2 at N 128 and 256 (so the 2D advantage of PURE at t = T comes from transport). Falsifier: ratio at most 1.0 at N 128 and 256.
- GATE (fixed before the runs): a candidate passes if its 1 - n.n_e L1 is at most 1.05 times plicRDF's on the 3D sphere at N 64 and 128 in both positions.

No transport run with a candidate is started; Tomislav decides that from the gate.

### Part B: bounds control (print only, existing keys)

deformation3D N 64, Co 0.5, serial, ringSumBoundTrace 1, clip false, snapTol 0, exactly the setup of job 57098598 (the same binary). Five runs: row (c) of the normal control (the PURE loop with normalModel plicRDF) with nAlphaBounds 3; plicRDF with nAlphaBounds 6 and 12; PURE with nAlphaBounds 6 and 12. Per run, next to the two existing runs at 3: min alpha and max alpha - 1 over the run before the bounding and after the last pass; the pass-by-pass minimum at the step of the lowest alpha; the five worst cells; the number of steps with min alpha < -1e-6; E_shape at the end; mass error; total wall clock.

- PB1: row (c) has plicRDF's bounds: min alpha before the bounding and after pass 3 within a factor 1.5 of plicRDF's (-0.0118 and -2.1e-4). Falsifier: after pass 3 below -1e-3.
- PB2: plicRDF converges with the passes: |min alpha| after the last pass at most 2e-5 with 6 passes and at most 1e-7 with 12. Falsifier: above 1e-4 with 6.
- PB3: PURE does not converge: |min alpha| after the last pass above 5e-4 with 12 passes. Falsifier: below 5e-5.
- PB4: E_shape of each row changes by less than 1 % between 3, 6 and 12 passes.

Tables: etc/normal_gate_tables.md and etc/bounding_pass_tables.md.

### Results of the static normal gate and the bounds control, 2026-11-23 (tables: `etc/normal_gate_tables.md`, `etc/bounding_pass_tables.md`; jobs: bounds 57105150, 57105152 to 57105155 (job 57098598 for the runs at 3 passes of plicRDF and PURE); gate 57105659 to 57105672 (the N 32 position 0 job also 57105546, without the plicRDF debug print), the position 1 jobs 57106311 to 57106317 (the first position 1 submissions, 57105660, 57105667 to 57105672 and 57106050 to 57106056, died in initDavofFields: the shifted centre was written with doubled parentheses by my sed; the table holds the reruns); build 57105297, commit dfbf795)

Implementation of the control model ringSumRDF (`src/VoF/davof/normalModels/ringSumRDFNormal.[HC]`): the ring sum of every interface cell, the plane positioned by `davof::positionPlane` (one `friend` declaration in `davof.H`), then rdfIterations times `reconstructedDistanceFunction::constructRDF` on the stencil of `markCellsNearSurf(interfaceCell, 1)` and the least-squares gradient (`leastSquareGrad polyDegree1`, the body of `plicRDF::gradSurf()`, which is private, over the davof scheme's interface cells), normalised, the planes repositioned between the repetitions, the scheme positioning the final plane. plicRDF's code is not touched: its `reconstruct()` starts from its own interpolated normals and its `gradSurf()` is private, so the model calls the public classes in plicRDF's order with plicRDF's arguments. The model's cells are the davof scheme's interface cells (face-fraction and ring-sum rule), plicRDF's are the cells with alpha inside (surfCellTol, 1 - surfCellTol); the cell counts agree in every state of the table. Citation: the header of plicRDF gives Scheufler and Roenby, J. Comput. Phys., 2019, doi 10.1016/j.jcp.2019.01.009 (no volume, pages or arXiv number); the header of ringSumRDF adds 383 (2019) 1-23 and arXiv:1801.05382, the same paper and DOI. The ring quadric exists on 826d8bb as the model `ringFit` (`ringFitNormal.[HC]`; the weighted least-squares Monge patch through the ring's plane centroids and unit normals; minCells 6, alignTol 0.2, weight area; ringArea normal as the fit fallback); it was run as it is, with every key explicit.

Part N (the static normal, one step of deltaT/1000 on the exact initial state; the 3D sphere at N 32/64/128, the 2D circle at N 64/128/256/512, each at the standard and the shifted centre; the shifted 2D circle is moved by 0.013 in x only, the y component of (0.013 0.007 0) lies along the circle's axis):
- rdfIterations 0 reproduces the ringSum row at N 32 (3D, standard centre): the last line of error.dat is identical in every column.
- 1 - n.n_e L1 of plicRDF, 3D sphere, standard centre N 32/64/128: 4.16e-4, 1.47e-4, 3.33e-5; shifted: 6.39e-4, 1.50e-4, 3.53e-5. Ratio of ringSum to plicRDF: 1.90, 2.06, 2.14 (standard) and 1.79, 1.93, 2.06 (shifted). Ratio of ringSumRDF with one iteration: 1.176, 1.046, 1.019 (standard), 1.069, 1.018, 1.024 (shifted); with two iterations: 1.131, 0.997, 0.977 and 1.035, 0.978, 0.983. 2D circle, ratio of ringSumRDF with one iteration at N 64/128/256/512: 1.109, 0.984, 1.088, 0.990 (standard), 1.016, 1.059, 1.025, 0.981 (shifted); with two iterations 0.96, 0.90, 0.97, 0.89 and 0.90, 0.95, 0.93, 0.89.
- The ring quadric (ringFit), 1 - n.n_e L1 against plicRDF's, 3D sphere N 32/64/128: standard 0.022, 0.006, 1.78; shifted 0.365, 1.73, 0.362. Its 1 - n.n_e Linf is 1.2e-4, 2.0e-5, 1.37e-2 (standard) and 4.0e-2, 0.24, 2.3e-2 (shifted); |n - n_e| L1 is 0.13, 0.076, 0.33 (standard) and 0.25, 0.23, 0.12 times plicRDF's. The fit fell back to the ringArea normal in 0, 0, 57 cells (standard) and 1, 3, 7 cells (shifted) of N 32/64/128; the large Linf values were not investigated. On the 2D circle every interface cell falls back (76 of 76 at N 64, 612 of 612 at N 512): the ringFit rows equal the ringSum rows to the digits of the table.
- plicRDF's own RDF loop in this step (DebugSwitches plicRDF 1, row (a); two reconstruct calls print, the value is the loop index at the exit, iterations used = index + 1): 1 and 1 (3D N 32, 64, both positions), 2 and 1 (3D N 128), 2 and 1 (every 2D N, both positions). The logs of the shear3D matrix runs do not print plicRDF's iteration count; no rerun.
- PG1 is false: the ring quadric's ratio is at most 0.5 in 2 of the 4 cases (N 64 standard 0.006, N 128 shifted 0.362) and above 1.0 in 2 (N 128 standard 1.78, N 64 shifted 1.73); the falsifier (above 1.0 in any of the four) is met.
- PG2 is false: ringSumRDF with one iteration is within 10 % of plicRDF's 1 - n.n_e L1 in 12 of the 14 state/position/N cases (all but 3D N 32 standard, 1.176, and 2D N 64 standard, 1.109); beyond 25 % in none, so the falsifier is not met.
- PG3 holds: the ratio of |n - n_e| L1 of the ring-sum normal to plicRDF's on the 2D circle is 1.73 (N 128) and 1.94 (N 256) at the standard centre, 1.94 and 1.83 at the shifted centre (at least 1.2 in 4 of 4; 1 - n.n_e L1 ratios 3.1 to 3.8 and 2.9 to 3.2).
- GATE (at most 1.05 times plicRDF's 1 - n.n_e L1 on the 3D sphere at N 64 and 128 in both positions): ringSum no (1.93 to 2.14); the ring quadric no (0.006, 1.78, 1.73, 0.362); ringSumRDF with one iteration yes (1.046, 1.019, 1.018, 1.024); with two iterations yes (0.997, 0.977, 0.978, 0.983). No transport run was started with a candidate.

Part B (the bounds control, deformation3D N 64, target Co 0.5, serial; solver Courant maximum 0.4950 in every run; the binary of job 57098598 in a copy of its platforms directory):

| row | passes | min alpha before the bounding | after the last pass | steps with min alpha < -1e-6 | E_shape (end) | wall (s) |
|---|---|---|---|---|---|---|
| plicRDF | 3 / 6 / 12 | -0.0118 / -0.0157 / -0.0157 | -2.1e-4 / -5.2e-8 / -7.0e-13 | 67 / 0 / 0 | 2.979e-3 / 3.019e-3 / 3.019e-3 | 321 / 316 / 320 |
| PURE | 3 / 6 / 12 | -0.0287 / -0.0297 / -0.0297 | -5.3e-3 / -2.3e-17 / -2.3e-17 | 265 / 0 / 0 | 3.860e-3 for all three | 287 / 287 / 285 |
| (c) PURE loop, plicRDF normal | 3 | -0.0139 | -2.7e-4 | 56 | 3.133e-3 | 335 |

- Pass by pass at the step of the lowest alpha: plicRDF with 6 passes (step 55) -0.00808, -1.74e-3, -4.6e-4, -4.0e-5, -1.9e-5, -4.8e-6, -5.2e-8; with 12 (step 897) -1.5e-3, -8.5e-8, -1.9e-11, -7.0e-13 (it stops at 4 passes); PURE with 6 and 12 passes (step 860) -0.0225, -1.0e-4, -8.3e-5, -3.5e-5, -2.3e-17 (it stops at 4 passes, the runs are identical); row (c) at 3 passes (step 59) -0.00644, -2.3e-3, -7.5e-4, -2.7e-4. The five worst cells of each run are in the table.
- PB1 holds: row (c) has min alpha -0.0139 before the bounding (1.18 times plicRDF's) and -2.7e-4 after pass 3 (1.32 times); not below -1e-3.
- PB2 holds: plicRDF's |min alpha| after the last pass is 5.2e-8 with 6 passes (at most 2e-5) and 7.0e-13 with 12 (at most 1e-7).
- PB3 is false: PURE's |min alpha| after the last pass is 2.3e-17 with 12 passes (and with 6); the falsifier (below 5e-5) is met. With 3 passes it is 5.3e-3.
- PB4 is false for plicRDF and holds for PURE: plicRDF's E_shape changes by +1.35 % from 3 to 6 passes and -0.010 % from 6 to 12; PURE's by -0.001 % and 0.000 %.
- Mass error (max |dev|) 1.5e-11 % to 3.7e-11 % in all runs.

## Branch davof-isoAdvector-normal, 2026-11-23, next batch: the transport matrix with HYB1 at 6 bounding passes and the ring quadric's outliers; pre-registration (Tomislav, written before the first run)

Rules: same workspace and branch, tmaric only, OpenFOAM Co 0.2 and 0.5 only with the solver's Courant line checked, every davof key explicit, no repair, no tuning, no new model.

### Part M: transport matrix

One binary for every row: the current build of the branch (with ringSumRDF), its hash recorded in the table. Cases: shear2D N 64/128/256 (ranks 1/1/8), shear3D and deformation3D N 32/64/128 (ranks 1/8/64), Co 0.2 and 0.5: the 18 cases of the curved matrix, same dictionaries. Rows, all with nAlphaBounds 6 (the only change against the curved matrix), clip false, snapTol 0: (1) plicRDF; (2) PURE (normalModel ringSum, key set as before); (3) HYB1 = PURE's key set with normalModel ringSumRDF, rdfIterations 1. One repetition per case, except the six N 128 cases on 64 ranks of the 3D flows, which get three repetitions per row (median wall clock reported, all three values). Columns as in entering_matrix_tables.md (corrected columns) plus: wall-clock ratio to plicRDF at equal ranks, E_shape ratio to plicRDF, and for rows (1) and (2) the E_shape ratio to the existing 3-pass value of the curved matrix.

- PH1: HYB1's E_shape is within 3 % of plicRDF's in six of six shear3D cases and within 6 % in six of six deformation3D cases. Falsifier: beyond 10 % in any 3D case.
- PH2: HYB1 saves little: its median wall clock at N 128 on 64 ranks is at least 0.85 times plicRDF's in four of four 3D cases. Falsifier: at most 0.70 in any of them.
- PH3: on shear2D HYB1 loses PURE's advantage: HYB1 / plicRDF in E_shape at least 0.9 at N 64, both Co. Falsifier: at most 0.75.
- PH4: with 6 passes every run of every row has min alpha >= -1e-6 and max alpha - 1 <= 1e-6. Falsifier: any PURE run below -1e-4.
- PH5: against the 3-pass values of the curved matrix, PURE's E_shape changes by less than 0.1 % and plicRDF's by less than 2 % in every 3D case.

### Part Q: the ring quadric's outliers (print only, static, no transport)

State: the 3D sphere, shifted centre, N 64 (ratio 1.73, Linf 0.24), and for comparison the standard centre N 64 (ratio 0.006). For the ten cells with the largest 1 - n.n_e of ringFit in each state: alpha_k, |m_k|/A_k, the number of ring planes that enter the fit, whether the fit fell back, whatever conditioning measure the fit already computes (no solver is added; if none exists, say so), 1 - n.n_e of ringFit and of ringSum in that cell. Also the count of interface cells with 1 - n.n_e above 1e-3 for ringFit, ringSum and plicRDF in both states, and ringFit's L1 with these ten cells removed.

- HQ1: the outliers are sliver cells: at least eight of the ten worst cells of the shifted state have alpha_k outside [0.01, 0.99]. Falsifier: fewer than five.
- HQ2: without the ten worst cells ringFit's 1 - n.n_e L1 in the shifted state is at most 0.1 times plicRDF's. Falsifier: above 0.5.

Tables: etc/transport_hyb_tables.md and etc/quadric_outlier_tables.md.

### Results of the transport matrix with HYB1 at 6 passes and of the ring quadric's outliers, 2026-11-23/24 (tables: `etc/transport_hyb_tables.md`, `etc/quadric_outlier_tables.md`; jobs: matrix 57106834 to 57106851 (all 18 completed, no abort), outlier study 57107251, 57107252, build of the print-only keys 57106972)

Part M. One binary for every row (the build at commit dfbf795, libVoF sha256 c1da5c86ddf849a9, copy `d2_iso6.sbatch` of the iso-time route, rows RDF, PURE, HYB1, nAlphaBounds 6, REPS 1, 3 for the six N 128 cases on 64 ranks); all 54 runs completed; solver Courant labels in the table. Output `iso_hyb/`.
- PH1 is false: HYB1 / plicRDF in E_shape is 1.047, 1.039, 1.030, 1.027, 1.028, 1.020 on shear3D (N 32/64/128, Co 0.2 then 0.5; within 3 % in 3 of 6 as the table counts, the 1.030 case at the limit) and 1.018, 1.007, 1.254, 1.273, 1.034, 1.010 on deformation3D (within 6 % in 4 of 6); beyond 10 % in 2 of 12 (deformation3D N 64, both Co: 1.254, 1.273), so the falsifier is met.
- PH2 is false and its falsifier is not met: the median wall clock of HYB1 over plicRDF at N 128 on 64 ranks is 0.816, 0.776 (shear3D Co 0.2, 0.5), 0.752, 0.763 (deformation3D): at least 0.85 in 0 of 4, at most 0.70 in 0 of 4. PURE over plicRDF is 0.597, 0.572, 0.527, 0.541. The three repetitions agree within 3 % except in deformation3D, where one value of HYB1 (745 s against 658, 661 at Co 0.2) and of plicRDF (934 against 874, 879) lies 6 to 13 % above the others.
- PH3 holds: HYB1 / plicRDF in E_shape on shear2D N 64 is 0.941 (Co 0.2) and 1.019 (Co 0.5).
- PH4 is false and its falsifier is not met: 33 of 54 runs have min alpha >= -1e-6 and max alpha - 1 <= 1e-6. Over all runs min alpha is -4.8e-5 (plicRDF), -3.4e-5 (PURE), -4.5e-5 (HYB1) and max alpha - 1 is 1.6e-4, 4.2e-6, 7.1e-5; no PURE run is below -1e-4. The exceptions are the N 128 cases of the 3D flows (64 ranks), shear2D N 128 Co 0.5 and N 256 (8 ranks) and deformation3D N 64 (8 ranks) PURE; the same case deformation3D N 64 Co 0.5 PURE with 6 passes ended at -2.3e-17 in the serial run of the bounds control and ends at -4.9e-6 in this 8-rank run; the cause was not investigated.
- PH5 holds: against the 3-pass values of the curved matrix PURE's E_shape changes by at most 0.006 % in all 12 3D cases; plicRDF's by at most 1.35 % (deformation3D N 64 Co 0.5; the others below 0.22 %) in all 12.

Part Q (the ring quadric's outliers, static, 3D sphere N 64; print-only keys `normalCellDiag` of reconstructionError and `ringFitDiag` of ringFit in a separate build of this branch, commit eb27227; the fit computes no conditioning measure, only the success of its Cholesky solve, and none was added):
- Shifted centre: 1 - n.n_e L1 plicRDF 1.50e-4, ringSum 2.88e-4, ringFit 2.60e-4 (1.73 times plicRDF's, as in the gate); cells above 1e-3: plicRDF 49, ringSum 143, ringFit 20; ringFit's fallbacks 3 cells. The ten worst ringFit cells have 1 - n.n_e 0.237, 0.050, 0.035, 0.032, 0.024, 0.012, 0.011, 0.0092, 0.0078, 0.0061; all have alpha_k within 2e-5 of 1 and |m_k|/A_k between 4e-6 and 2.5e-4, 13 to 15 ring planes entering the fit (9 to 15 with n.e3 > alignTol), none fell back; ringSum's and plicRDF's errors in the same cells are 1e-4 to 2.4e-3. ringFit's L1 without these ten cells: 1.58e-5, 0.105 times plicRDF's.
- Standard centre: ringFit's L1 9.2e-7 (0.006 times plicRDF's); cells above 1e-3: plicRDF 0, ringSum 216, ringFit 0; the ten worst ringFit cells (all alpha_k 0.99984, 1 - n.n_e at most 2.0e-5); no fallback.
- HQ1 holds: 10 of 10 of the ten worst cells of the shifted state have alpha_k outside [0.01, 0.99] (all nearly full cells with a tiny interface, |m_k|/A_k at most 2.5e-4).
- HQ2 is false by 5 %: ringFit's L1 without the ten worst cells is 0.105 times plicRDF's (the limit 0.1); the falsifier (above 0.5) is not met.

## Branch davof-isoAdvector-normal, 2026-11-24, closing items: the table correction, the ring quadric under transport, the lower bound in parallel; pre-registration (Tomislav, written before the first run of Parts R and P)

Rules unchanged: same workspace, tmaric only, OpenFOAM Co 0.2 and 0.5 only, own job ids only, every davof key explicit, no repair, no tuning, no new model.

### Part C: correct transport_hyb_tables.md

Its normal columns are shifted again (the column "abs(n-n_e) L1" reads the mean cell volume, the column "1-n.n_e Linf" carries the corrected tables' 1 - n.n_e L1): regenerate with the corrected indices, check three rows against entering_matrix_tables.md, name the wrong script in the README, and check that normal_gate_tables.md and quadric_outlier_tables.md do not share the defect.

### Part R: the ring quadric under transport (existing model, no code change)

Row QUAD = PURE's key set with normalModel ringFit exactly as run in the gate, nAlphaBounds 6. Cases: shear3D and deformation3D, N 32/64/128 (ranks 1/8/64), Co 0.2 and 0.5: 12 runs, one repetition, on the Part M binary (libVoF c1da5c86...). Reported: the Part M columns (corrected), the E_shape and wall-clock ratios to plicRDF and to PURE, the fit fallbacks per step (maximum and at the last step) against the interface cell count.

- PR1: all 12 runs complete with min alpha >= -1e-4 and max alpha - 1 <= 1e-4. Falsifier: an abort, or min alpha below -1e-3 in any run.
- PR2: the quadric does not beat plicRDF under transport, because the transported face fractions and not the fit limit the normal: QUAD / plicRDF in E_shape on shear3D at N 128 lies between 1.0 and 1.4 at both Co. Falsifier: below 0.95 at both Co.
- PR3: fit fallbacks at the last step are below 5 % of the interface cells in all six shear3D runs.
- Decision rule, fixed now: if PR2's falsifier is met, the quadric is a lead; otherwise the normal question is closed in the fork.

### Part P: the lower bound in parallel (print only, existing trace)

deformation3D N 64, Co 0.5, PURE, nAlphaBounds 6 ended at -2.3e-17 in the serial bounds control and at -4.9e-6 on 8 ranks in Part M. Rerun on 8 ranks with ringSumBoundTrace 1 (if the trace does not work in parallel: say what blocks it and stop this part). Reported: min alpha before the bounding and after each pass at the step of the lowest final alpha; for the cells below -1e-9 after the last pass over the run: how many have a face on a processor patch, how many are within one cell of one, how many neither; the five worst cells with these flags.

- HP1: at least 80 % of the cells below -1e-9 after the last pass have a face on a processor patch or share a face with such a cell. Falsifier: fewer than 40 %.

Tables: etc/quadric_transport_tables.md and etc/bounding_parallel_tables.md.

### Results of the closing items, 2026-11-24 (tables: `etc/transport_hyb_tables.md` regenerated, `etc/quadric_transport_tables.md`, `etc/bounding_parallel_tables.md`; jobs: QUAD 57115466 to 57115477, trace in parallel 57115478)

Part C (the table correction). The script of the previous batch, `d2_hyb-tab.py`, read the reconstructionError columns with the uncorrected indices 11, 12, 13 (avgVol, LNormalDiff1, LNormalDiff2) instead of 12, 13, 14 (LNormalDiff1, LNormalDiff2, LNormalDiffInf): the column "abs(n-n_e) L1" showed the mean cell volume and the "1-n.n_e Linf" column carried 1 - n.n_e L1. It was written from a stale local copy of `d2_ent-tab.py`, the one from before that script's own column correction. The table is regenerated with the indices 12, 13, 14; three rows checked against `entering_matrix_tables.md`: shear3D N 128 Co 0.2 plicRDF (E_shape 3.307e-4, |n - n_e| L1 0.053, 1 - n.n_e L1 0.00232, Linf 0.0157 in both) and PURE (5.035e-4, 0.0688, 0.00358, 0.0251 in both), and deformation3D N 64 Co 0.5 PURE (|n - n_e| L1 0.711 against 0.709 and 1 - n.n_e L1 0.432 against 0.431 with 6 against 3 passes; E_shape 3.860e-3 in both). The other tables do not share the defect: `normal_gate_tables.md` is built by `d2_gate-tab.py` with the indices 12, 13, 14 (the gate's plicRDF 1 - n.n_e L1 at N 64, standard centre, 1.47212e-4, equals LNormal1 of error.dat, 1.47212e-4; the ring quadric, shifted centre, 2.59508e-4 equals LNormal1 too), and `quadric_outlier_tables.md` reads the per-cell lines of the print-only key (its totals agree with the gate's to all digits shown). The E_shape-based statements PH1 to PH5 and their numbers are unchanged (E_shape does not come from these columns).

Part R (the ring quadric under transport; row QUAD = PURE's key set with normalModel ringFit, applyTo all, weight area, alignTol 0.2, minCells 6, nAlphaBounds 6; shear3D and deformation3D, N 32/64/128, Co 0.2 and 0.5; the Part M binary, libVoF c1da5c86ddf849a9; one repetition). All 12 runs completed.
- E_shape of QUAD / plicRDF: shear3D 1.264, 1.355, 1.454 (Co 0.2) and 1.270, 1.357, 1.431 (Co 0.5) at N 32/64/128; deformation3D 1.144, 1.011, 1.221 and 1.181, 1.047, 1.210. QUAD / PURE: shear3D 0.947, 0.968, 0.955 and 0.936, 0.963, 0.953; deformation3D 0.910, 0.788, 1.106 and 0.959, 0.819, 1.107. Wall clock of QUAD / plicRDF: shear3D 2.73, 1.28, 0.887 and 1.94, 0.908, 1.01; deformation3D 2.10, 1.30, 0.846 and 1.60, 1.19, 0.854 (QUAD / PURE 1.49 to 3.11).
- PR1 is false and its falsifier is met: 8 of 12 runs have min alpha >= -1e-4 and max alpha - 1 <= 1e-4; exceptions shear3D N 64 Co 0.2 (min -2.0e-4), shear3D N 128 Co 0.5 (min -7.1e-3, max alpha - 1 1.9e-3), deformation3D N 64 Co 0.5 (-1.5e-4), deformation3D N 128 Co 0.5 (max alpha - 1 6.8e-4); min alpha below -1e-3 in 1 run (shear3D N 128 Co 0.5). That run has 133591 interface cells at the last step.
- PR2 is false and its falsifier is not met: QUAD / plicRDF in E_shape on shear3D at N 128 is 1.454 (Co 0.2) and 1.431 (Co 0.5), just above 1.4 at both; below 0.95 at neither. Decision rule: the falsifier is not met, so the normal question is closed in the fork.
- PR3 is false: the fit fell back to the ringArea normal in 1612 of 2354 (N 32 Co 0.2), 26056 of 31190 (N 64), 3893 of 11694 (N 128) and 3598 of 4535, 2783 of 4884, 113486 of 133591 (Co 0.5) interface cells at the last step of the six shear3D runs (33 % to 85 %); below 5 % in 0 of 6. The maximum over the steps is in the table.

Part P (the lower bound in parallel; deformation3D N 64, Co 0.5, PURE, nAlphaBounds 6, 8 ranks, `ringSumBoundTrace 1`). The trace runs in parallel but only the pass-by-pass minima are global; the counts, the cells and the five worst cells it prints after the loop are those of the master rank (local labels, no reduction, no processor-patch information; in 11 of the 960 steps the master holds no alpha in the printed range), so the per-cell flags the part asks for cannot be produced without changing the trace. Stopped there. Global pass by pass at the step of the lowest final alpha (step 768): before -1.9e-3; after passes 1 to 6 -7.6e-4, -7.6e-4, -3.6e-5, -2.1e-5, -1.0e-5, -4.9e-6; the bounding is not converged after the 6 passes in 21 of the 960 steps, where the serial run converges to round-off in at most 4. HP1 not evaluated.
