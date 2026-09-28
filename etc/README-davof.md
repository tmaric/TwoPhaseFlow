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

## How leia finds it

leia's `etc/leia-env.sh` appends `$LEIA_TPF_DIR/platforms/<WM_OPTIONS>/{bin,lib}` after
leia's own directories; the default is the `TwoPhaseFlow` clone next to the leia clone.
The `interFlow` benchmark arm of leia's studies (`config/interFlowDroplet{2D,3D}.yaml`)
therefore runs this build without any study-config reference to this repository.
