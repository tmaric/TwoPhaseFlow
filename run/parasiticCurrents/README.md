# interFlow aborts at exit with a corrupted heap (OpenFOAM-v2512)

**Status:** reproducible, cause unknown, blocks using `interFlow` results as
publication-grade numbers. Two minimal cases in this directory reproduce it in
seconds.

---

## The symptom

Every `interFlow` run completes its time loop normally, prints `End`, and *then*
aborts during static destruction. Exit code **134** (SIGABRT). glibc reports one of

```
malloc_consolidate(): invalid chunk size
malloc_consolidate(): unaligned fastbin chunk detected
corrupted size vs. prev_size
```

The stack at the abort is in library teardown, not in the solver:

```
__libc_free
Foam::dictionary::~dictionary          (libOpenFOAM.so)
__cxa_finalize
Foam::debug::deleteControlDictPtr::~deleteControlDictPtr
```

So the heap was already corrupted while the solver ran; glibc only detects it when
it consolidates the heap at exit. **The results are written, but nothing certifies
them** — an out-of-bounds write happened somewhere, and where it landed is unknown.

## Which case to run

```bash
source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc
# plus whatever exports put this repository's build on PATH/LD_LIBRARY_PATH

cd run/parasiticCurrents/stationaryDroplet2D
./Allrun
```

Takes a few seconds (64x64x1 cells, 20 steps). `Allrun` runs `interFlow` directly
rather than through `runApplication`, because the exit code is the whole point and
`runApplication` swallows it. Expected output:

```
interFlow exit code: 134
steps completed:     20
reached 'End':       1
glibc heap message:  corrupted size vs. prev_size
```

`stationaryDroplet3D` is the same physics in 3D (60^3, 20 steps, ~1 min) and behaves
identically. Start with 2D.

## What has been ruled out

| Hypothesis | Test | Result |
|---|---|---|
| MPI teardown artefact | ran serial and at np=4 | **both** abort |
| A specific curvature model | swept all five `surfaceTensionForceModel`s: `RDF`, `fitParaboloid`, `gradAlpha`, `heightFunction`, `constantCurvature` | **all five** abort identically |
| Something in our case setup | ran the repo's own `run/velocityStaticCircle` with its own `blockMesh` + `initAlphaField` + `interFlow`, nothing changed but a shortened `endTime` | **also aborts**: 20 steps, `End`, `corrupted size vs. prev_size`, exit 134 |
| Wrong OpenFOAM version | `.github/workflows/openfoam.yml` pins `OF_VERSION: "2512"`; README states master compiles v2406–v2512 | v2512 is a **targeted** version |
| A build problem | `Allwmake` exit 0, no errors, `ldd` on `interFlow` has no missing entries, `interFlow -help` runs | build is clean |

So it is **not** in the curvature models, **not** MPI, **not** our cases, and **not**
a version mismatch. It is in the common path — the solver, `isoAdvection`/`plicRDF`,
or the `surfaceForces`/`deltaFunction` infrastructure — and it is pre-existing.

Environment where this was observed: OpenFOAM-v2512 (`Build: _87ed40d256-20251219`),
gcc 11.5.0, twophaseflow at `de9826f9ffb24f4b635ac97fd388ebd560cfc174`
("Merge branch 'pr-61'"). Reproduced on both a WSL/Ubuntu laptop and the
Lichtenberg cluster (Red Hat, gcc 11.5.0, OpenMPI 4.1.8), i.e. it is not
machine-specific.

## Suggested next step

Rebuild with AddressSanitizer and run the 2D case — that should name the offending
write directly:

```bash
export WM_COMPILE_OPTION=Debug     # or add the flags to your wmake rules
# add to c++FLAGS / linker flags:  -fsanitize=address -fno-omit-frame-pointer
./Allwmake
cd run/parasiticCurrents/stationaryDroplet2D && ./Allrun
```

ASAN aborts at the *first* invalid access with a full stack, rather than at exit
like glibc does. Note the solver is heavily templated OpenFOAM code, so an ASAN
build of the whole repo plus `libOpenFOAM` may be needed if the report points into
OpenFOAM itself. `valgrind --tool=memcheck` on the 2D case is a slower alternative
that needs no rebuild.

Two smaller things worth checking first, both cheap:

- `heightFunction.C:103` has an unused variable `avgColVal`; the file does raw
  stencil indexing and is a plausible place for an off-by-one — though note the
  corruption also occurs with `constantCurvature`, which never enters it.
- The repo builds its own `libVoF` alongside OpenFOAM's `geometricVoF`. If any
  class is defined in both with different layouts, a cross-library allocation and
  free would corrupt the heap exactly like this. Worth checking for duplicate
  symbol names between `libVoF.so` and OpenFOAM's `libgeometricVoF.so`.

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
at the same point. That is exactly why the heap corruption needs resolving rather
than working around.

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

One initialisation trap worth recording: a 2D droplet must be initialised as
`type cylinder`, not `type sphere`. On a one-cell-thick mesh `sphere` integrates a
sphere *slab*, giving `sum(alpha*V) = 9.5260e-10` against the correct
`pi R^2 t = 9.8175e-10` — 3% low — and the wrong Laplace jump (`2 sigma/R` instead
of `sigma/R`). `initAlphaField` accepts `plane`, `sphere`, `cylinder`, `paraboloid`,
`sin`, `composedFunction` and `ellipsoidImplicitFunction` (note that last name).
