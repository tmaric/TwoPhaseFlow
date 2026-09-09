# Acoustic methods-paper submission configuration

The submission uses **geometric PLIC face averaging with the established
one-field acoustic flux**. It does not use the experimental multipoint
transmission flux.

Put the following in `system/fvSchemes`:

```foam
acousticInterface
{
    areaFraction       plicAverage;
    flux               legacy;
    writeDiagnostics   false;
}
```

The library and solver need to be built once. The word `legacy` above selects
the retained flux; the area calculation is geometric. The original method is
still available by setting both entries to `legacy`. Missing entries retain
their original defaults for compatibility.

## Case mapping

The following files explicitly select the submission method:

| Case | Configuration source |
|---|---|
| Homogeneous plane waves | `run/acousticTests/FrequencyDomainTests/homogeneousPlaneWave2D/system/fvSchemes.in`; case preparation generates `fvSchemes` from this template |
| Layered gas/liquid interface and PML | `run/acousticTests/FrequencyDomainTests/layeredInterface1D/system/fvSchemes` |
| Baffled-piston radiation | `run/acousticTests/FrequencyDomainTests/pistonRadiation/system/fvSchemes` |
| Resolved rigid sphere / Gorkov force | `run/acousticTests/FrequencyDomainTests/Andrade2019/WedgeLevitatorDropGorkov/system/fvSchemes` |
| Fixed reconstructed droplets | `run/acousticTests/FrequencyDomainTests/Andrade2019/WedgeLevitatorDropAlpha/system/fvSchemes` |

The fixed-droplet case supports the associated software/application workflow.
Other solvers and moving-interface cases are not changed by this selection.

## Method and supported claims

For each internal or processor face, use both valid neighbouring PLIC cuts
when available, averaging their unsigned liquid fractions; use a single valid
cut directly. A zero/full cut remains a valid contribution. The calculation
does not use `phi` to select a donor. Canonical triangulation makes the
nonplanar-face calculation consistent under orientation reversal. Same-phase
pure cells use zero/one; opposite pure cells without a plane retain linear
interpolation. Missing mixed-cell planes and invalid geometry are errors.

The pressure flux coefficient remains
`a_f = 1/(alpha_f*rho_l + (1-alpha_f)*rho_g)`.
Cell density, reaction quadrature, PML coefficients, boundary conditions,
nonorthogonal correction, and velocity/force postprocessing are unchanged.
Liquid acoustic transmission remains part of the continuum model. This
coefficient closure is not an exact sharp-interface reconstruction of
phase-wise pressure gradients.

Controlled comparisons support a reproducibility benefit, not a general
pressure-accuracy improvement:

| Finest controlled comparison | Original areas + original flux | Geometric areas + original flux |
|---|---:|---:|
| Planar transmission sequences | Baseline errors | Identical within roundoff |
| 48³ penetrable-sphere relative complex-pressure error | 0.02974595 | 0.02975350 |
| Same sphere: serial/eight-rank relative pressure difference | 2.6395e-6 | 4.9163e-15 |
| Same sphere: stored matrix entries | 3,041,280 | 3,041,280 |
| Finest fixed spherical-droplet reflector force | 1.32623e-5 | 1.32623e-5 |
| Finest deformed-droplet reflector force | 5.49337e-6 | 5.42176e-6 |

These auxiliary sphere and droplet comparisons are documented in
[the verification report](testsuite/acousticInterface/VERIFICATION.md).
Their pressure/force differences should not be described as evidence of
improved accuracy. The existing limitations for shifted interfaces and
strongly nonorthogonal meshes remain.

## Reproduction and verification

Build instructions are in [INSTALL_ACOUSTICS.md](INSTALL_ACOUSTICS.md).
The working branch is `feature/acousticLevitation`. The implementation
was developed from revision `5b1d46437356dc4debfdd8e33932d201e3c63c88`;
that base revision alone does not contain the implementation changes.
Use the accompanying committed implementation and its source manifest.

For the current WSL checkout, select the isolated tested binaries with:

```sh
source .review-testing/interface-implementation/env.sh
```

Build the compiled check utility after the solver, and use a Python environment
with NumPy and pytest installed (see INSTALL_ACOUSTICS.md):

```sh
wmake apps/benchmark/testAcousticInterface
```

Run the compiled/integration tests:

```sh
python3 -m pytest --confcutdir=testsuite/acousticInterface \
    testsuite/acousticInterface/test_invariants.py
```

Check that the selected paper dictionaries preserve the existing layered-PML
results in fresh directories, against the frozen original executable:

```sh
python3 testsuite/acousticInterface/submission_check.py \
    --output /absolute/new/submission-pml-comparison \
    --baseline /absolute/path/acousticHelmholtzFoam-baseline \
    --ranks 1,4
```

This exercises both 8,000-cell material orderings, the five aligned mesh levels
560/1120/2240/4480/8960, and weak PML damping. It inherits the explicit
submission dictionary, prepares each case once, compares complex fields to
the original executable, checks finite output and assembled residuals, and
retains all raw cases and manifests. General three-method comparisons use a
neutral prepared dictionary and apply each selection afterwards, so changing
the submission templates does not redefine the baseline.

The manuscript update is made in `helmHoltzPML-revised.tex` in the article
workspace. Original manuscript sources and published figure/table inputs are
preserved. The revised text defines the candidate cuts, geometric averaging,
pure-cell fallback, processor exchange and retained flux explicitly. Internal
review notes remain in the source behind a disabled conditional and are
excluded from the submission PDF.

Implementation source and original benchmark evidence are indexed in
`testsuite/acousticInterface/results/final-provenance.json` and
`testsuite/acousticInterface/results/verification-summary.json`.
The configuration-update check is recorded separately in
`testsuite/acousticInterface/results/submission-selection.json`.

## Verification of the selected configuration (9 September 2026)

All 20 compiled/integration pytest tests passed in WSL and on Lichtenberg.
The layered-PML check completed 32 serial/four-rank comparisons locally and
32 serial/eight-rank comparisons in cluster job **54506190**. The maximum
pressure difference from the original executable was 1.85e-12 locally and
1.95e-12 on the cluster, below the 1e-9 gate.

Reapplying the paper's original 1,200-point VTK sampling to the selected-method
fields reproduced its reported values: 1.1970121e-3 for the 8,000-cell forward
case, 4.4245418e-3 in reverse, and 2.6568582e-4 at the finest aligned resolution.
All probe flags were valid. The postprocessor now fails on missing/invalid
samples or nonfinite values rather than silently discarding them. An outside-
domain probe was explicitly checked to fail.

Direct cell-centre errors differ from these sampled errors (for example,
7.8441e-5 at N=8960); the manuscript now specifies the sampling operation so
these two metrics are not confused. Existing plotted curves and table values
are preserved. Raw cases, logs and manifests are retained under
`.review-testing/submission-geometric-averaging` on both hosts.

## Completed manuscript evidence audit

Section 4.4 and Table 4 now contain the quantitative comparison of the two area definitions. All 19 figures and 7 tables have a recorded source. The audit corrected stale homogeneous and piston assets and regenerated layered assets from the selected-method fields. See [the manuscript provenance report](testsuite/acousticInterface/MANUSCRIPT_PROVENANCE.md) for the complete mapping, 12 additional damping runs, and the limits of retained historical force records.
