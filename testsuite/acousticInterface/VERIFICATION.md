# Selectable acoustic interfaces: implementation and verification

**Submission decision (9 September 2026):** use `plicAverage/legacy`.
See [../../SUBMISSION_ACOUSTICS.md](../../SUBMISSION_ACOUSTICS.md).
The failed experimental transmission-flux gates below do not describe a
replacement selected for the submission; that prototype remains opt-in.


The three requested combinations are implemented and selectable at startup,
using the same solver and library binaries. Missing configuration preserves
the legacy method. **The full scientific acceptance gate has not passed.**
The new flux is suitable for controlled experiments, but the present evidence
does not support recommending it as the default or claiming general second-order
accuracy.

The most consequential results are:

1. Forward, high-contrast transmission with a quarter-cell interface offset
   becomes substantially less accurate with the new flux, despite a correctly
   assembled and accurately solved linear system.
2. Both the retained nonorthogonal iteration and the broader skewed-mesh
   accuracy gate fail on some oblique meshes. Failed runs are retained.
3. The new flux substantially improves penetrable-sphere pressure errors and
   reversed layered transmission.
4. Geometric averaging removes decomposition dependence observed in the
   preserved legacy area method. The original behavior is intentionally still
   selectable.

Implementation reference date: 2026-09-09. Raw results, fields, prepared meshes,
logs and failed cases are under `.review-testing/interface-implementation/`
in the WSL and Lichtenberg checkouts. The compact result table and figure are
in [results/verification-summary.json](results/verification-summary.json) and
[results/pressure-comparison.svg](results/pressure-comparison.svg).

## Acceptance status

| Gate | Status and evidence |
|---|---|
| Build and selection | Passed on OpenFOAM v2606 in WSL and v2512 on Lichtenberg. All three pairings run without rebuilding. Unknown methods, invalid numerical controls and `legacy/plicTransmission` are rejected. |
| Preserve legacy | Passed against an independently rebuilt original executable. Missing dictionary and explicit legacy settings reproduce the baseline within `1e-9`; the serial comparisons agree to roundoff. The two legacy area/exchange headers remain byte-identical to the originals. |
| Compiled local geometry and flux | Passed 4,036 checks, including 780 affine-fit configurations. Maximum affine-flux absolute error on these scaled test configurations: `2.9421e-15`. The configured bounds are `1e-12` for controlled geometry and `1e-10` for local flux checks. |
| Dictionary, I/O and integration tests | 20 pytest tests passed both locally and on the final cluster build, including homogeneous liquid and gas with nonzero PML coefficients. |
| Actual 1D assembled operator | Passed comparison with an independently assembled path-resistance matrix, constant preservation, expected Neumann nullspace, Dirichlet nullspace removal and a checkerboard test. |
| Actual 3D assembled operator | Passed on 512 cells: direct face-operator integration versus PETSc matrix action `3.55e-16`; unscaled conservation and constant-state defects `1.53e-15`. |
| Small diffusion spectrum | Passed on the tested 1D and 3D matrices. The 3D Neumann operator has one numerical null mode and no negative-real mode outside roundoff; the Dirichlet operator has positive real eigenvalues. Its relative symmetry defect is `0.0189`, so symmetry/positive definiteness is not assumed. |
| New-method MPI consistency | Passed on completed controlled cases at local ranks 1/2/4 and cluster ranks 1/2/4/8. The perturbed oblique study's maximum accepted transmission difference was `1.91e-12`; the full sphere study's was `3.01e-15`. |
| Shared processor operator | The 48³ sphere audits cover 241 processor faces locally (stencils reaching three ranks) and 386 faces on eight cluster ranks (stencils reaching four ranks), with exactly opposite exported weights. The solver also checks processor scalar areas against `1e-12`. |
| Legacy MPI consistency across all meshes | Failed on some sphere/oblique cases, reproducing an existing limitation. At 48³ the eight-rank legacy pressure difference is `2.64e-6`; geometric averaging reduces it to `4.92e-15`. The maximum legacy difference across the sphere sequence is `2.97e-4`. |
| Linear-solve evidence | Completed verification cases require positive PETSc convergence reason, finite solution, true relative assembled residual at most `1e-10`, and reported backward error. Failed residual gates remain failed. |
| Nonorthogonal iteration | Failed on selected oblique cases. A small linear residual does not imply convergence of the outer correction. Drivers require the final relative pressure change at most `1e-8` when multiple corrections are requested. |
| Finite-frequency accuracy | Mixed; see results below. The forward quarter-shift case rules out a blanket acceptance claim for the new flux. |
| General second-order claim | Not established. Aligned planar cases and reversed quarter-shift transmission reach approximately second order. The 3D sphere's measured new-method order is `1.69` and the axisymmetric sphere's is `1.74`. |
| Derived velocity and force accuracy | Not established by the pressure tests. Existing definitions are unchanged; selected diagnostics and force comparisons are reported separately. |

These are measured checks on the listed configurations, not a proof of global
stability for every admissible mesh or density contrast.

## Controlled implementation

See [../../src/acousticInterface/README.md](../../src/acousticInterface/README.md)
for configuration, limitations and build instructions.

- `src/acousticInterface` contains startup selection, legacy wrapping,
  geometric cuts, bounded graph-halo exchange and cached face operators.
- `transmissionFit.C` contains canonical triangulation and triangle clipping,
  scaled exact owner/neighbour constraints, a nullspace/SVD least-squares fit,
  and elimination into global-cell pressure weights.
- The solver masks the old physical Laplacian coefficient on replaced faces,
  removing both its implicit and explicit nonorthogonal contributions there.
  It inserts the replacement real weights in both pressure diagonal blocks.
- Each unscaled interior face contributes equal and opposite rows, followed by
  the existing cell-density scaling. Processor faces have a single operator
  builder; the opposite operator is exchanged.
- PML coupling, reaction quadrature, density interpolation, boundary
  definitions and velocity/force postprocessing are retained.
- Additional local and remote columns are included in PETSc preallocation.
  The allocation estimate is a safe superset of the actual row stencil;
  `MAT_NEW_NONZERO_ALLOCATION_ERR` remains enabled.

PETSc distinguishes locally owned and remote columns in parallel preallocation.
Its convergence reason is positive for convergence and negative for divergence.
The implementation uses these interfaces and additionally evaluates the actual
residual. [PETSc preallocation](https://petsc.org/release/manualpages/Mat/MatMPIAIJSetPreallocation/),
[PETSc convergence reason](https://petsc.org/release/manualpages/KSP/KSPGetConvergedReason/).

The C++ changes after the large refinement studies tightened rejection of
processor-cyclic and overset patches and added cell/face/plane details to error
messages. They did not change supported-case numerical weights or assembly.
The boundary-guard source passed baseline/3D MPI gates in job `54503962`;
job `54504147` rebuilt the final diagnostics source and passed the full test
suite and baseline/eight-rank regression. Source manifests distinguish these
builds from the earlier study binaries. The final compiled source hashes match
between hosts; see [results/final-provenance.json](results/final-provenance.json).

## Layered transmission

The domain is `[0,0.14] m`, frequency 20 kHz, with analytical complex Dirichlet
pressure at both ends and no PML. The interface is at
`0.07 + offset*(0.14/N)`. Material properties are `rho=1.2, c=343` and
`rho=1000, c=1500` in SI units; the reverse study swaps their ordering.

Both studies contain 108 runs: six sizes `160,240,320,480,640,960`,
three offsets `0,0.25,0.5`, three methods and ranks 1/8.
Uniform cell volumes make the reported norm identical to the volume-weighted
complex-pressure relative L2 norm. Phase and near-interface errors are also
retained. Here, “forward” means gas-to-liquid ordering.

| N | Legacy / geometric areas, forward offset 0.25 | Transmission flux, forward offset 0.25 |
|---:|---:|---:|
| 160 | 0.152380 | 1.65304 |
| 240 | 0.0616428 | 2.32027 |
| 320 | 0.0343425 | 4.96009 |
| 480 | 0.0161444 | 2.28082 |
| 640 | 0.00997201 | 0.753569 |
| 960 | 0.00548084 | 0.259184 |

The new-method error is nonmonotone on the coarse sequence. A fit using only
the finest three points gives slope 3.10, but this is not evidence of an
asymptotic third-order method: it follows a large erroneous response.

An independent 1D matrix using exact geometric path resistances reproduces
the computed pressure to relative `7.18e-14` in the audited failure case.
This rules out PETSc insertion as the explanation for that particular error.
The interaction of an affine transmission flux with finite-frequency phase
curvature and the retained mixed-cell reaction quadrature is a plausible
explanation, not a demonstrated complete diagnosis. No unrequested quadrature
change was made to improve these results.

| Finest N=960 case | Legacy relative error | New relative error | New fitted order, finest three |
|---|---:|---:|---:|
| Forward, aligned | 0.00298567 | 0.00298567 | 2.018 |
| Forward, half-cell offset | 0.00354558 | 0.00354558 | 2.070 |
| Reverse, aligned | 0.000102101 | 0.000102101 | 2.001 |
| Reverse, quarter-cell offset | 0.00770812 | 0.000102981 | 2.012 |
| Reverse, half-cell offset | 0.000103668 | 0.000103668 | 2.021 |

For reverse quarter-shift transmission, the legacy fitted order is 1.039.
The geometry-only method agrees with legacy in this 1D setting.

`metrics-supplement.json` files add reflection/transmission magnitudes and
phases, and velocity errors, from the original retained fields without
rerunning or overwriting the studies. Phases are in radians; the “gas” and
“liquid” field labels follow the configured phase names, whose material values
are swapped in the reverse experiment.

The forward `N=960` quarter-shift case has velocity relative L2 errors
`0.00253` with legacy and `0.10779` with transmission. Pure-phase errors
also differ substantially. These diagnostics use the original gradient-based
velocity calculation, not the new face flux as a replacement velocity.

## Penetrable sphere

The independent reference evaluates a regular interior spherical-Bessel
expansion and an outgoing exterior spherical-Hankel expansion, matched by
pressure and inverse-density normal-derivative continuity. It uses 25 angular
orders, `a=0.001 m` and exterior `ka=1`, with the same liquid/gas properties
as above. A separate numerical check verifies the reference's interface
conditions.

The full 3D domain is the cube `[-2a,2a]^3`, with analytical outer pressure
data. Meshes are `16³,24³,32³,48³`; each was run with all three methods,
serial/four ranks in WSL and serial/eight ranks on Lichtenberg.

| N | Legacy error | Geometric areas + legacy flux | Transmission flux |
|---:|---:|---:|---:|
| 16 | 0.921181 | 0.911995 | 0.00874945 |
| 24 | 0.179596 | 0.179828 | 0.00388465 |
| 32 | 0.0647163 | 0.0647667 | 0.00254319 |
| 48 | 0.0297459 | 0.0297535 | 0.00120921 |

The new method improves the finest-grid error by about a factor of 25.
Its fitted order over `24³,32³,48³` is 1.694. Phase-weighted and interface
errors are retained, including the finest new liquid error `0.004643` and
gas error `0.0005947`.

The axisymmetric reference uses a one-degree wedge, radial extent `2a` and
axial extent `[-2a,2a]`, with `N` radial and `2N` axial cells. Sizes
`16,32,64,96` were run locally with all methods at ranks 1/4. At `N=96`,
legacy/geometric errors are approximately `0.003635` and transmission error
is `0.00010617`. The transmission order over the finest three is 1.744.
This verifies the centre-plane reconstruction with actual wedge face areas.

### Direct-solver cost

For the 48³ cluster case on eight ranks:

| Method | Stored matrix nonzeros | MUMPS summed factor memory, MB | Comparison wall time, seconds |
|---|---:|---:|---:|
| Legacy | 3,041,280 | 3,680 | 23.86 |
| Geometric areas + legacy flux | 3,041,280 | 3,525 | 23.23 |
| Transmission flux | 3,181,584 | 3,881 | 28.06 |

Memory is MUMPS `INFOG(22)`, not total process peak memory. Timings include
the driver's decomposition/solve/reconstruction work and are single
measurements, not a performance scaling study. The new stencil increases
stored nonzeros by about 4.6% and factor memory relative to the geometric-only
case by about 10%. The direct-solver bottleneck remains.

## Oblique interfaces and iteration failures

The oblique tests use a fixed physical plane and 2D meshes sheared relative to
that plane by 15, 45 and 75 degrees. This keeps the interface away from
nonsymmetry physical boundaries while retaining exact analytical data.
Two mesh variants are provided: uniform shear (“regular”) and a sinusoidally
varying shear (“perturbed”). The latter fixes the end boundaries.
The transverse mesh has eight cells; longitudinal sizes are `48,96,192`.
These are directional refinement studies, not an isotropic 2D convergence
proof.

At 2 kHz with 100 nonorthogonal corrections, the perturbed cluster study
contains 54 runs. Twenty-two failed: all 75-degree runs and the coarsest
45-degree legacy/geometric runs. The successful new-method error sequences
are:

| Angle | N=48 | N=96 | N=192 |
|---|---:|---:|---:|
| 15 degrees | 0.00394327 | 0.00137540 | 0.000308218 |
| 45 degrees | 0.00728295 | 0.00145377 | 0.000269944 |

The local regular-shear study also contains 54 runs. All 45- and 75-degree
runs failed the outer-iteration gate (36 failures); all 15-degree runs
completed. The new 15-degree errors are `0.00515583,0.00205438,0.000360243`.
The matching eight-rank cluster study produced the same 36 failures out of
54 runs. Its maximum accepted new-method MPI difference was `2.66e-13`.

The original 20 kHz oblique pilot failed to converge and is preserved.
Reducing the frequency for the controlled comparison does not erase that
failure. These results do not support a claim that the new local flux fixes
the solver's nonorthogonal iteration on arbitrary skewed meshes.

The homogeneous 2 kHz wave comparison has 18 completed runs (three sizes,
three methods, ranks 1/8), with maximum serial/parallel difference
`9.32e-13`. Separate local homogeneous liquid/gas tests exercise nonzero
PML coupling and equality across the three selections. These PML checks
establish compatibility, not a new PML accuracy/reflection study.

## Existing wedge droplets and forces

The comparison driver prepares new copies of the existing
`Andrade2019/WedgeLevitatorDropAlpha` study using cfMesh and wedge extrusion.
Nominal mesh sizes are `0.0004` and `0.0002 m`; sphere and aspect-ratio-two
droplets use the same prepared state for each method. The meshes contain
3,355 and 13,190 cells.

There are 24 successful local comparisons (ranks 1/4) and 24 successful
cluster comparisons (ranks 1/8). The retained coded reflector force output,
in the existing wedge convention, is:

| Shape and mesh size | Legacy | Geometric areas | Transmission flux |
|---|---:|---:|---:|
| Sphere, 0.0004 | 1.20884e-5 | 1.20382e-5 | 1.39642e-5 |
| Sphere, 0.0002 | 1.32623e-5 | 1.32623e-5 | 1.43655e-5 |
| Aspect 2, 0.0004 | 4.50280e-6 | 4.47417e-6 | 8.01843e-6 |
| Aspect 2, 0.0002 | 5.49337e-6 | 5.42176e-6 | 7.74529e-6 |

No independent exact force is available for these application comparisons.
The force change is not evidence of improved accuracy. Pressure/velocity RMS,
raw fields and `pr`/`momFlux` finite checks are retained.

## Reproduction and provenance

WSL repository: `/home/tmaric/OpenFOAM/repos/TwoPhaseFlowAcoustics`.
Cluster repository: `/work/scratch/tm83tomy/TwoPhaseFlowAcoustics`.
The studies ran on `feature/acousticLevitation`, starting from
`5b1d46437356dc4debfdd8e33932d201e3c63c88`, with the implementation
identified by the recorded source hashes. This is the historical build
identity; the implementation is now retained in this branch's commit history. The subsequent methods-paper update and its data audit are documented in [MANUSCRIPT_PROVENANCE.md](MANUSCRIPT_PROVENANCE.md).

The original WSL executable SHA-256 is
`6d701dd6925966f14ed12b9afed1aa9687abdcd6914bd3e82f704e0ce2a8ed20`.
The reference was built from copied original sources before introducing the
new library. Builds use isolated `build/bin` and `build/lib` directories;
normal shared OpenFOAM executables are not replaced.

Comparison drivers require fresh output directories. They snapshot inputs
after preparation, apply dictionary selections to separate run copies, record
source and binary hashes, and check that binaries are unchanged throughout
each comparison. The final driver additionally archives the relevant source
files. Earlier study snapshots and transfer archives remain available beside
their manifests. Missing fields, NaNs, empty data and incomplete runs fail.
There are no probe-validity flags in these comparisons because raw cell fields
are used directly.

| Cluster job | Work and recorded outcome |
|---|---|
| 54503491 | Initial environment/build setup failed before testing. |
| 54503505 | Dependent study cancelled after that setup failure; no study results. |
| 54503619 | Corrected build, compiled checks and baseline/3D gates completed. |
| 54503635_0 | Forward and reversed layered studies completed, 216 total runs. |
| 54503635_1 | Perturbed oblique and homogeneous studies produced 72 rows; 22 oblique failures. This earlier driver returned success despite failed rows; the current driver returns failure when any row fails. |
| 54503635_2 | Full 3D sphere sequence completed, 24 runs. |
| 54503636, 54503712 | Wedge setup failed because the cluster cfMesh plugin directory contained no sources. |
| 54503713 | Wedge build and 24 comparisons completed after supplying the WSL cfMesh source snapshot. |
| 54503962 | Final boundary-guard source rebuilt; 18 pytest tests and baseline/sphere gates on ranks 1/2/4/8 completed. |
| 54504068 | Regular-shear oblique study completed all 54 attempts; 36 numerical failures are retained and propagated to the job exit status. |
| 54504147 | Final diagnostics source rebuilt; all 20 pytest tests, baseline/eight-rank regression and the retained 48³ processor-operator audit passed. |

Cluster builds used GCC 11.5, OpenMPI 4.1.8, OpenFOAM v2512 and PETSc 3.25.
Local builds used OpenFOAM v2606. The final supported-case numerical source
is identical on both hosts; host-specific binary hashes are recorded separately.

Useful verification commands, after sourcing the documented environment:

```sh
python3 -m pytest --confcutdir=testsuite/acousticInterface \
    testsuite/acousticInterface/test_invariants.py

python3 testsuite/acousticInterface/audit_operator.py \
    --prepared /absolute/prepared-sphere-N8 \
    --output /absolute/new-operator-audit

python3 testsuite/acousticInterface/audit_parallel.py \
    --case /absolute/parallel-transmission-case

python3 testsuite/acousticInterface/advanced.py \
    --output /absolute/new-oblique-study --scenario oblique \
    --mesh-style regular --sizes 48,96,192 --angles 15,45,75 --ranks 1,4

python3 testsuite/acousticInterface/summarize.py \
    --results .review-testing/interface-implementation \
    --output testsuite/acousticInterface/results
```

## Limits requiring further work before a general recommendation

The finite-frequency shifted-interface failure needs a separate investigation
of phase curvature, reaction quadrature and interface location. Any proposed
change should remain separately selectable so this comparison remains
interpretable. The nonorthogonal iteration also needs its own stability work.

The exact-plane tests in this implementation exercise the compiled local
geometry/fit functions; global cases use the existing PLIC reconstruction.
An independent prescribed-plane override for global PDE studies has not been
added. The test set does not exhaust every concave face, reconstruction
setting, mesh degeneracy or contact geometry. Axisymmetric and Cartesian
sphere refinements are bounded, and the oblique refinements are directional.

The new modes reject unsupported coupled/overset patches, and transmission
rejects a changing mesh, interface/PML overlap and mixed-cell support at
nonsymmetry physical boundaries. Rank-deficient or ill-conditioned fits stop
with diagnostics. General contact-boundary treatment, AMI/cyclic support,
moving interfaces and a scalable segregated Helmholtz solution remain outside
this implementation. `ROADMAP.md` remains locally ignored.
