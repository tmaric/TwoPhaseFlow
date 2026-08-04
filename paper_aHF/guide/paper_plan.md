# Preliminary Split: Software Paper vs Numerical Method Paper

This note divides the current `acousticHelmholtzFoam.tex` material according to
the annotations in `ahf_commented.pdf`. The main decision is whether the current
manuscript should be framed as a CPC-style software paper, a numerical-method
paper, or split into two related manuscripts.

## High-Level Decision

The current manuscript should be split conceptually into:

1. A software paper centered on the OpenFOAM implementation, usability,
   verification cases, performance, and reproducibility.
2. A numerical method paper centered on the governing equations,
   heterogeneous/two-phase material modeling, finite-volume discretization,
   interface coefficient treatment, PML formulation, consistency, and accuracy.

The same repository and cases can support both papers, but the emphasis and
section order should differ.

## Paper A: Software / CPC Paper

### Possible Title

`acousticHelmholtzFoam`: An OpenFOAM module for parallel frequency-domain
acoustic simulations in heterogeneous media

Alternative if we name the broader project:

`TwoPhaseFlow-acoustics`: OpenFOAM tools for frequency-domain acoustics in
heterogeneous two-phase configurations

### Main Claim

The paper presents an open-source OpenFOAM implementation for solving
frequency-domain acoustic Helmholtz problems in heterogeneous media, with PML
domain truncation, MPI/PETSc block-system assembly, verification cases, example
applications, and performance measurements.

### Keep / Emphasize

- `Program summary`.
- Repository, license, dependencies, build instructions, and tested environment.
- Solver overview: fields, dictionaries, workflow, inputs, outputs.
- `setPMLFields` utility and supported PML configurations.
- OpenFOAM implementation details:
  - real/imaginary pressure storage;
  - volume and surface coefficient fields;
  - block-coupled PETSc matrix assembly;
  - processor-patch handling;
  - output fields such as `pa`, `Ure`, `Uim`, `pr`, and `momFlux`.
- Verification cases as reproducible tutorial/regression cases.
- Application example: wedge levitator/drop case.
- Performance/scaling section.
- Practical limitations and extension points.

### Shorten / Move Out

- Full derivation of the heterogeneous wave equation.
- Long algebraic derivation of the PML coefficients.
- Detailed numerical-method novelty discussion.
- Deep discussion of interface coefficient consistency, except as a brief
  implementation note and reference to the method paper.

### Proposed Section Order

1. Introduction
2. Software Overview
3. Mathematical and Numerical Background
   - short governing equation summary only;
   - cite/point to method paper for full derivation if available.
4. Implementation in OpenFOAM
   - solver fields and dictionaries;
   - material-property construction;
   - PML utility;
   - PETSc/MPI assembly.
5. Verification and Example Cases
   - piston radiation;
   - one-dimensional layered interface;
   - wedge levitator/drop example.
6. Performance
7. Availability, Build, and Reproducibility
8. Conclusions

### Current Manuscript Mapping

- `Program summary`: software paper.
- `Introduction`: reuse, but frame as software contribution.
- `Mathematical model`: compress heavily.
- `Numerical implementation/Solver overview`: software paper.
- `Numerical implementation/PML configuration`: software paper.
- `Verification`: software paper, presented as reproducible tests.
- `Application example`: software paper.
- `Performance`: software paper.
- `Conclusions`: rewrite around software contribution and repository.

## Paper B: Numerical Method Paper

### Possible Title

A finite-volume method for frequency-domain acoustic wave propagation in
heterogeneous two-phase media with interface-aware material coefficients and PML
truncation

### Main Claim

The paper presents and analyzes the finite-volume formulation for the
heterogeneous Helmholtz equation with two-phase material coefficients,
interface-aware face coefficients, real/imaginary block coupling, and PML terms.
The contribution must be stated as a numerical method, not only an
implementation.

### Keep / Emphasize

- Governing equations and perturbation derivation.
- Frequency-domain Helmholtz formulation for heterogeneous media.
- Acoustic radiation-pressure derivation, if it is part of the method/result.
- PML coordinate stretching and coefficient derivation.
- Two-phase material coefficient modeling:
  - linear VOF density;
  - linear acoustic compressibility `1/(rho*c^2)`;
  - relation between `rho`, compressibility, and `k^2`;
  - face coefficient treatment for fluxes.
- Finite-volume discretization:
  - cell-centered operator;
  - face interpolation;
  - non-orthogonal correction;
  - boundary-gradient treatment.
- Accuracy/consistency questions:
  - what order is expected;
  - how boundary gradients are evaluated;
  - whether observed convergence is supported by enough mesh levels.
- PML sensitivity and reflection/attenuation behavior.
- Numerical verification focused on method accuracy, not only software tests.

### Shorten / Move Out

- Program summary.
- Build instructions and dependency details.
- PETSc/MUMPS branding unless solver choice affects the numerical result.
- Detailed OpenFOAM dictionary descriptions.
- Performance/scaling, unless relevant to algorithmic complexity.
- Software architecture and repository organization.

### Proposed Section Order

1. Introduction
   - state numerical-method gap and novelty.
2. Governing Equations and Frequency-Domain Model
3. Two-Phase Material Representation
   - one-field density and compressibility;
   - interface-cell and face coefficient treatment.
4. PML Formulation for the Heterogeneous Helmholtz Equation
5. Finite-Volume Discretization
   - diagonal and coupling operators;
   - real/imaginary block form;
   - boundary and non-orthogonal gradient treatment.
6. Verification
   - analytical piston case;
   - layered-interface case;
   - mesh-convergence study with enough resolutions.
7. Discussion
   - limitations;
   - boundary-gradient accuracy;
   - what has and has not been proven/tested.
8. Conclusions

### Current Manuscript Mapping

- `Introduction`: reuse, but reframe around numerical gap/novelty.
- `Mathematical model`: numerical method paper.
- `Acoustic radiation pressure`: numerical method paper if used as a derived
  quantity; otherwise shorten.
- `Perfectly matched layer`: numerical method paper.
- `Finite-volume discretisation`: numerical method paper, but update the
  density/compressibility model.
- `Solver overview`: mostly software paper; keep only what is needed to define
  the algebraic system.
- `PML configuration`: mostly software paper; keep only mathematical profile
  definitions if needed.
- `Verification`: numerical method paper, but strengthen convergence claims.
- `Performance`: software paper, not numerical paper.

## Immediate Manuscript Actions

1. Update the current two-phase coefficient equations to use the corrected
   linear VOF density and linear acoustic compressibility.
2. Decide whether the current target remains CPC. If yes, keep the current file
   as the software paper and move the longer numerical derivations into a
   separate method-paper draft or appendix.
3. If a numerical paper is planned, create a second manuscript file with the
   method-focused outline above.
4. Rewrite the abstract according to the chosen paper type:
   - software paper: emphasize OpenFOAM module, reproducible cases, PML utility,
     PETSc/MPI implementation, and availability;
   - numerical paper: emphasize method, discretization, coefficient treatment,
     PML formulation, and accuracy.
5. Revisit convergence and boundary-gradient claims before making strong
   numerical-method statements.

| Case | Numerical paper | Software paper |
|---|---|---|
| Homogeneous plane-wave convergence | Full orthogonal/warped convergence study | One compact regression result only |
| Layered-interface convergence | Full interface-position study | Remove detailed derivation and convergence |
| Isolated PML reflection | Full parameter and resolution study | Remove sensitivity study |
| Integrated interface + PML | Full analytical comparison | Keep as one end-to-end software example |
| Baffled piston | Full near/far-field verification | Remove |
| `setPMLFields` profile check | Not needed | Add |
| Serial/MPI solution equivalence | Not needed | Add |
| Levitator height sweep | Not needed | Keep as application |
| Strong/problem-size scaling | Not needed | Keep |