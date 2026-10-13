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

Nothing changes for the other reconstruction schemes (plicRDF after ten
steps of the deformation case is bit-identical with the build before this
change, alpha and interfaceNormal). Dictionary keys of the alpha solver:
`storeAreaFractions` (default: the scheme is davof), the davof scheme's
`wispTol` (default 1e-3 since this step, leia's value), `isoFaceTol`.

### The mixed normal and the normal-source controls (2026-10-13)

The closed loop (own normal from alphaf, the alpha-positioned plane, alphaf
from the moved plane) amplifies through sliver interface cells (small
|m_c| compared with h^2: the positioned plane's sensitivity h/A, the own
normal's 1/A). Following Basilisk's mixed Youngs-centred scheme (Aulisa et
al. 2007, `src/myc.h`) the `davof` scheme keeps the slivers as interface
cells and gives them a stencil normal; keys of the alpha solver dictionary
(`davof.H`):

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

## The advection benchmarks (run/benchmark/advection)

- `planeTranslation3D`: the exact-solution gate of the face update (above).
- `sphereTranslation3D`: a sphere (R = 0.15) carried by U = (1 0 0) over
  one unit (3.3 diameters) in the 2 x 1 x 1 box at Courant 0.5, N = 32, 64,
  128 cells across the height; rows isoAdvection + davof on the
  `initDavofFields` state (quadraticFaces, detrixheAslam), + plicRDF and
  + isoAlpha on the setAlphaField state; `areaFractionError` every step,
  the reconstruction error at the end, `volumeFractionError` every step.
  README there for the measured numbers: with the rules above the davof
  loop (own normal from the transported face fractions, the plane
  positioned by alpha, the plane's fluxes) is UNSTABLE on the sphere at
  N = 32 and 64: the face-fraction error grows from 1.2e-2 after one step
  to 0.4 within ten, alpha leaves [0, 1] at step 4, and the interface-cell
  count explodes; plicRDF and isoAlpha on the same runs are healthy. The
  deformation regression (`testsuite/advection/test-deformationSphere` at
  32 x 32 x 64, davof quadraticFaces, serial) is unbounded as well
  (Ebound 1.4 at t = 3, Emass 4e-3, Eshape 3.8e-2; plicRDF: 2e-8, 3e-16,
  4.5e-3). The gate attributes the growth to the sliver interface cells
  (the C2 mechanism); which rule of the loop changes is Tomislav's
  decision (the proposal's Stage 1 projection, or a plane-setting threshold
  above the interface-cell threshold).

## How leia finds it

leia's `etc/leia-env.sh` appends `$LEIA_TPF_DIR/platforms/<WM_OPTIONS>/{bin,lib}` after
leia's own directories; the default is the `TwoPhaseFlow` clone next to the leia clone.
The `interFlow` benchmark arm of leia's studies (`config/interFlowDroplet{2D,3D}.yaml`)
therefore runs this build without any study-config reference to this repository.
