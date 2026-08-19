# interFlow aborts at exit: duplicate statics in libVoF and libgeometricVoF (OpenFOAM-v2512)

**Status: FIXED.** Root cause was a **double `free()` of one static `Foam::word`
during library teardown**, caused by TwoPhaseFlow's `libVoF` and OpenFOAM-v2512's
`libgeometricVoF` both defining the same static objects. The forked classes now live
in an inline namespace `Foam::twoPhaseFlow`, so no symbol collides; see
[Fixing it](#fixing-it). `interFlow` exits 0 and is AddressSanitizer-clean, and the
solution is bit-identical to before the fix.

**The results are valid.** ASAN reports exactly one memory error in a full run, and
it happens *after* `End`, during static destruction. There are **zero** buffer
overflows or use-after-free events during the time loop, and the solution is
bit-identical to a run that does not abort. The earlier suspicion that "the heap was
already corrupted while the solver ran" was wrong.

---

## The symptom

Every `interFlow` run completes its time loop, prints `End`, and *then* aborts during
static destruction with exit code **134** (SIGABRT), glibc reporting one of

```
corrupted size vs. prev_size
malloc_consolidate(): invalid chunk size
malloc_consolidate(): unaligned fastbin chunk detected
```

The gdb stack is in library teardown, not in the solver:

```
__libc_free  <-  primitiveEntry::~primitiveEntry  <-  dictionary::~dictionary
             <-  debug::deleteControlDictPtr::~deleteControlDictPtr  <-  __cxa_finalize
```

## Root cause

`libVoF.so` (this repository) and `libgeometricVoF.so` (OpenFOAM-v2512) both define
**the same 19 static data symbols**, because this repository carries a fork of
OpenFOAM's `geometricVoF`:

```
Foam::reconstructionSchemes::typeName          Foam::reconstructionSchemes::debug
Foam::reconstruction::plicRDF::typeName        Foam::reconstruction::plicRDF::debug
Foam::reconstruction::isoAlpha::typeName       Foam::reconstruction::isoAlpha::debug
Foam::reconstruction::gradAlpha::typeName      Foam::reconstruction::gradAlpha::debug
Foam::cutCell::debug                           Foam::cutFace::debug
Foam::reconstructionSchemes::componentsConstructorTablePtr_        ... and 8 more
```

The dynamic linker collapses each pair to **one address**, but **each library keeps
its own static constructor and its own `__cxa_finalize` destructor registration** for
it. So the object is constructed twice and, at exit, **destroyed twice** — `~word()`
frees the same heap pointer twice. That corrupts glibc's free lists, and glibc
notices later while consolidating the debug-switch dictionary.

ASAN names it exactly:

```
ERROR: AddressSanitizer: attempting double-free on 0x5030001a2c90
  #1 __cxa_finalize
  #2 libVoF.so+0x6f1f6                          <-- second free, from libVoF
0x5030001a2c90 is located 0 bytes inside of 22-byte region
freed by thread T0 here:      __run_exit_handlers                <-- first free
previously allocated by thread T0 here:
  #1 Foam::word::word(char const*, bool)  in libgeometricVoF.so  <-- allocation
```

The 22-byte region is `"reconstructionSchemes"` — 21 characters plus NUL — i.e.
`Foam::reconstructionSchemes::typeName`.

### Why `libgeometricVoF` is in the process at all

`interFlow` does not link it. It arrives through the function objects:

```
controlDict: libs (fieldFunctionObjects)
  -> libfieldFunctionObjects.so  ->  libreactingMultiphaseSystem.so
  -> ... -> libincompressibleMultiphaseSystems.so  ->  libgeometricVoF.so
```

### Why the abort looked erratic

The double-free happens **whenever both libraries are loaded**, but glibc only
*aborts* when the freed chunk happens to sit next to one it later validates. Function
objects change the heap layout, so they change whether glibc notices — they do not
cause the bug. Confirmed: a run with `-noFunctionObjects` plus `libs (geometricVoF)`
exits 0 **and still contains the double-free** under ASAN.

## Evidence

| Experiment | Result |
|---|---|
| Baseline, function objects on | exit 134, `corrupted size vs. prev_size` |
| `-noFunctionObjects` (`libgeometricVoF` never loads) | **exit 0**, ASAN clean |
| `-noFunctionObjects` + `libs (geometricVoF)` | exit 0 but **ASAN still reports the double-free** |
| `LD_PRELOAD=libgeometricVoF.so` (reverses symbol binding) | exit **139**, segfault — same bug, different manifestation |
| Only `fieldMinMax`; only `volFieldValue` | each aborts alone |
| `probes` from `libsampling` (already linked, no `geometricVoF`) | exit 0 |
| Merely `dlopen`ing `fieldFunctionObjects` with no FO running | exit 0 |
| `libVoF` relinked `-Wl,-Bsymbolic` (DF_SYMBOLIC verified set) | **still exit 134** — cannot help, each library's destructor list runs regardless |
| ASAN, full run | **1** error, at teardown; **0** overflows/UAF during the time loop |
| Solver output, aborting run vs clean run | **identical** (same md5 over all solver lines) |

Previously ruled out and still ruled out: MPI (serial aborts too), any specific
`surfaceTensionForceModel` (all five abort), our cases (`run/velocityStaticCircle`
aborts too), version mismatch, and a broken build.

## A second, latent hazard in the same fork

Independent of the abort, `plicRDF` has a **different layout** in the two libraries:

```
OpenFOAM-v2512   reconstructedDistanceFunction  RDF_;   sizeof(plicRDF) = 2784
this repository  reconstructedDistanceFunction& RDF_;   sizeof(plicRDF) = 1936
```

848 bytes apart (856-byte object vs 8-byte reference), and `sIterPLIC_` sits after
it, so its offset differs too. Both libraries register `plicRDF` into the *same*
interposed `componentsConstructorTablePtr_`; whichever registers first wins, and that
depends on load order. Today `libVoF` loads first and stays self-consistent, but if
that order ever flips, a real in-run heap overflow becomes possible. This is not the
cause of the current abort, but it should not be left standing.

## Fixing it

The fork is real, not stale — 17 classes in `src/VoF` are duplicated and **all** have
diverged (`reconstructedDistanceFunction.C` by 428 lines, `isoAdvection.C` by 303,
`plicRDF.C` by 183), and the `RDF_` reference is a deliberate design difference. So
the two cheap options were both closed off:

- **Dropping the fork** and using OpenFOAM's `geometricVoF` is not a drop-in; the
  implementations and the `plicRDF` API have diverged.
- **Hiding the symbols** (`-fvisibility=hidden`, version script) breaks the build:
  `libsurfaceForces`, `libphaseChange` and `libpostProcess` all reference
  `Foam::reconstructionSchemes::typeName` as an undefined symbol.
- **`-Bsymbolic` / `-Bsymbolic-functions`** was tested and is ineffective (table
  above): it changes symbol *binding*, but both libraries still register a
  destructor for the same object.

### What was done

The forked classes are wrapped in an **inline namespace** `Foam::twoPhaseFlow`:

```cpp
namespace Foam
{
inline namespace twoPhaseFlow
{
namespace reconstruction
{
class plicRDF : public reconstructionSchemes { ... };
}
} // End inline namespace twoPhaseFlow
} // End namespace Foam
```

An inline namespace changes the **mangled** names — `plicRDF::typeName` becomes
`Foam::twoPhaseFlow::reconstruction::plicRDF::typeName`, which no longer collides —
while leaving every *spelling* of the name valid in the enclosing namespace. So
`Foam::reconstructionSchemes`, `reconstruction::plicRDF` and
`mesh.lookupObjectRef<reconstructionSchemes>(...)` keep compiling untouched
everywhere they already appear, in this repository and in downstream code. Nothing
outside the wrapped files had to change, and no `typeName` *string* changed, so case
dictionaries still say `reconstructionScheme plicRDF`.

Three families were wrapped, chosen by measuring which symbols actually collide
(strong `T`/`D`/`B` definitions shared with any OpenFOAM library) rather than by
guessing:

| Library | Forked classes wrapped | Collisions before → after |
|---|---|---|
| `libVoF` | `reconstructionSchemes`, `plicRDF`, `isoAlpha`, `gradAlpha`, `isoSurface`, `reconstructedDistanceFunction`, the `cutCell`/`cutFace` families, `surfaceIterator{Iso,PLIC}` | 97 → **0** |
| `libpostProcess` | `sampledInterface` | 20 → **0** |
| `libtwoPhaseModelThermo` | `phaseModel`, `twoPhaseMixtureThermo` | 55 → **0** |

`isoAdvection` needed nothing: it already sits in `Foam::advection` while OpenFOAM's
is `Foam::isoAdvection`, which is the same solution applied earlier by hand.

### Case-syntax change: `interface` → `reconstructedInterface`

Namespacing `sampledInterface` exposed a **second, independent conflict** that the
merged symbols had been masking. Both this repository and OpenFOAM's `geometricVoF`
register a `sampledSurface` under the *runtime name* `interface`, into OpenFOAM's
single `sampledSurface` table. That table keeps whichever library registered first
and rejects the second with

```
Duplicate entry interface in runtime table sampledSurface
```

Once the two classes were genuinely distinct, OpenFOAM won that key (it is loaded
first, indirectly via `libsampling`/`libfieldFunctionObjects`), so `type interface`
built *OpenFOAM's* sampler, which then looked for *OpenFOAM's* `reconstructionSchemes`
in the registry and found this repository's instead:

```
--> FOAM FATAL ERROR
    bad lookup of reconstructionScheme (objectRegistry fluid)
    expected a reconstructionSchemes, found a isoAlpha
```

This repository's sampler therefore registers under the unambiguous name
**`reconstructedInterface`**, and the 12 cases in `run/` and `testsuite/` that used
`type interface` were updated. **Existing external cases must make the same change:**

```
// before                    // after
type    interface;           type    reconstructedInterface;
```

Keeping `interface` as an alias was tried and dropped: OpenFOAM wins that key
whenever `libgeometricVoF` is loaded, so the alias never made the key usable and only
added a ten-line "Duplicate entry" stack trace to stderr on every run.

### Verification

| Check | Result |
|---|---|
| `Allwmake` from clean | exit 0, no errors; nothing outside the wrapped files changed |
| Strong-symbol collisions, every TwoPhaseFlow library vs **all** OpenFOAM libraries | **0** (except the turbulence library below) |
| `stationaryDroplet{2D,3D}`, `velocityStaticCircle`, `sinWave`, `stefanProblem`, `suckingInterface`, `cht/fixedFlux` | all exit **0**, all reach `End`, no `FOAM FATAL`, no `Duplicate entry`, no glibc message |
| ASAN (preloaded allocator) on all six cases | **0** errors — the double free is gone |
| Solver output, 2D case, before vs after the fix | **identical** (same md5 over all solver lines) |
| `max(p_rgh)` in the 2D case | 73.68 Pa against the exact Laplace jump 72.74 Pa |

### One known collision left

`libVoFphaseModelCompressibleTurbulenceModels` still shares 60 strong symbols with
OpenFOAM's `libVoFphaseCompressibleTurbulenceModels`. A namespace **cannot** fix
these: they are not forked classes but duplicate *explicit instantiations* of
OpenFOAM's own templates — `Foam::LESModel<...>`, `Foam::RASModels::kOmegaSST<...>`,
`Foam::laminarModels::Stokes<...>` and friends, instantiated for
`fluidThermoPhaseCompressibleTurbulenceModel` exactly as OpenFOAM's library does.
Moving them into `Foam::twoPhaseFlow` is not possible, because the templates belong
to OpenFOAM.

It is left alone because it is **not reachable**: walking the `DT_NEEDED` graph over
every OpenFOAM library, nothing depends on `libVoFphaseCompressibleTurbulenceModels`,
and in particular no function-object library pulls it in. Contrast
`libgeometricVoF`, which 19 libraries reach, three of them function-object libraries
(`libfieldFunctionObjects`, `liblagrangianFunctionObjects`, `libphaseFunctionObjects`)
— that reachability is exactly why this bug fired on ordinary cases. The turbulence
library can only be co-loaded by naming it explicitly in a `libs (...)` entry. The
real fix, if it ever matters, is to link OpenFOAM's library instead of
re-instantiating those templates.

## Reproducing

```bash
source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc
cd run/parasiticCurrents/stationaryDroplet2D
./Allrun
```

A few seconds (64x64x1 cells, 20 steps). `Allrun` runs `interFlow` directly rather
than through `runApplication`, because `runApplication` swallows the exit code.
To see the double-free itself:

```bash
ASAN_OPTIONS=detect_leaks=0:halt_on_error=0 \
  LD_PRELOAD=$(gcc -print-file-name=libasan.so) interFlow
```

No instrumented rebuild is needed — preloading ASAN's allocator is enough, because
both the allocation and the double free happen in library code.

`stationaryDroplet3D` is the same physics in 3D (60^3, 20 steps, ~1 min) and behaves
identically. Start with 2D.

Environment: OpenFOAM-v2512 (`Build: _87ed40d256-20251219`), gcc 11.5.0 (ASAN via
system gcc 13), twophaseflow at `de9826f` ("Merge branch 'pr-61'"). Originally
observed on both a WSL/Ubuntu laptop and the Lichtenberg cluster (Red Hat, gcc
11.5.0, OpenMPI 4.1.8), i.e. not machine-specific.

## Why we care

These cases are the stationary-droplet parasitic-current benchmark: a 1 mm droplet
of water in air, sigma = 0.07274 N/m, zero gravity, in a closed box. The exact
solution is **U = 0** with a piecewise-constant pressure, so any velocity that
appears is numerical. The exact Laplace jump is `sigma/R = 72.74 Pa` in 2D
(cylinder) and `2 sigma/R = 145.48 Pa` in 3D (sphere).

`interFlow` is being compared against OpenFOAM's `interFoam` and against a
semi-Lagrangian level-set solver, to separate *where the curvature comes from* from
*how the force is coupled*:

| solver | curvature source |
|---|---|
| `interFoam` | `kappa = -div(nHat)`, `nHat` from `interpolate(grad(alpha))` |
| level set | per-cell quadratic fit of a level set `psi` |
| **`interFlow` RDF** | signed distance **rebuilt from the PLIC planes** |
| **`interFlow` fitParaboloid** | paraboloid fit to the **reconstructed interface** |

All four are balanced-force CSF with the same `rAUf` cancellation, so the comparison
isolates the curvature source. `interFlow` matters most because it takes curvature
from a reconstructed geometric object — the same premise as the level-set method —
while conserving volume by construction. In the 20-step 2D case above it already
looks good: volume conserved to 12 digits, and `max|U|` about half of `interFoam`'s
at the same point.

## The full studies

The complete parameter studies live in the `leia` repository as snakemake configs,
not here, because they share that repository's study machinery:

- `config/interFlowDroplet2D.yaml` — N = 64/128/256, `STF_MODEL` in {RDF, fitParaboloid}
- `config/interFlowDroplet3D.yaml` — N_L = 60/76/95 on a 6 mm box, same two models
- `cases/interFlowDroplet{2D,3D}` — the case templates

Those are initialised with `leiaSetFields` instead of `initAlphaField`, so that all
three solvers start from a **bit-identical** volume fraction. The cases in this
directory deliberately avoid that dependency: they use this repository's own
`initAlphaField` so they can be run with nothing but OpenFOAM and this repo.

Note for those studies: `interFlow` now exits 0, so exit-code checks are meaningful
again. If those case templates sample the interface, they need
`type reconstructedInterface` rather than `type interface`.

One initialisation trap worth recording: a 2D droplet must be initialised as
`type cylinder`, not `type sphere`. On a one-cell-thick mesh `sphere` integrates a
sphere *slab*, giving `sum(alpha*V) = 9.5260e-10` against the correct
`pi R^2 t = 9.8175e-10` — 3% low — and the wrong Laplace jump (`2 sigma/R` instead
of `sigma/R`). `initAlphaField` accepts `plane`, `sphere`, `cylinder`, `paraboloid`,
`sin`, `composedFunction` and `ellipsoidImplicitFunction` (note that last name).
