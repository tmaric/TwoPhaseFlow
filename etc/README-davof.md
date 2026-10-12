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
   cell's plane sweeps it and sets `alphaf` AND `dVf` (the upwind value
   phi alpha_W dt equals the side value only for a pure upwind cell W, not
   for a sliver below the interface-cell threshold).
4. Boundary faces of a surface cell: the owner's plane sets `alphaf`; the
   outflow flux is the plane's as before, an inflow face keeps the boundary
   condition's flux.
5. The second ring: every face no plane has set whose cell shares a face
   with a surface cell is swept by the neighbouring plane whose
   end-of-step position is closest to the face centre (one and the same
   plane for a planar interface), which sets `alphaf` and `dVf` (inflow
   boundary faces excepted). A face the moving plane crosses during the
   step without belonging to a surface cell exists whenever a face-normal
   component of n is below the Courant number (the plane passes within
   Co h of a vertex ahead of both cells); without this rule the plane
   gate fails at the first step.
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

## How leia finds it

leia's `etc/leia-env.sh` appends `$LEIA_TPF_DIR/platforms/<WM_OPTIONS>/{bin,lib}` after
leia's own directories; the default is the `TwoPhaseFlow` clone next to the leia clone.
The `interFlow` benchmark arm of leia's studies (`config/interFlowDroplet{2D,3D}.yaml`)
therefore runs this build without any study-config reference to this repository.
