# Acoustic interface methods

`acousticHelmholtzFoam` selects the area and flux calculations at startup from
the optional `acousticInterface` dictionary in `system/fvSchemes`.

| `areaFraction` | `flux` | Meaning |
|---|---|---|
| `legacy` | `legacy` | Original donor/upwind PLIC areas and mixture-coefficient Laplacian |
| `plicAverage` | `legacy` | **Submission method:** symmetric geometric areas with the original flux |
| `plicAverage` | `plicTransmission` | Geometric areas with the experimental transmission flux |

Both entries default to `legacy`. Required inputs, including `alphaf`, remain compatible. The original area calculation and processor
exchange are preserved verbatim under `legacy/`. The solver's original copies
are retained for comparison. The combination `legacy/plicTransmission` and
unknown method names are rejected.

```foam
acousticInterface
{
    areaFraction       plicAverage;
    flux               legacy;
    writeDiagnostics   false;
}
```

Copyable snippets for all three combinations are under `examples/`.

Compile the library and solver once. Switching any supported combination then
requires only changing the dictionary and starting another run. Changing the
selection during an existing run is not supported.

## Build and use

Source the OpenFOAM environment and the repository's `scripts/bashrc`, then run:

```sh
bash testsuite/acousticInterface/build.sh
```

The root `Allwmake` also builds and links the new library. PETSc remains a
solver dependency; the acoustic-interface library has no PETSc dependency.
The implementation has been built against OpenFOAM v2606 in WSL and v2512 on
Lichtenberg.

For the isolated binaries produced during implementation in the WSL checkout:

```sh
source .review-testing/interface-implementation/env.sh
```

This selects the tested solver and library without replacing shared OpenFOAM
binaries. The original executable is preserved in
`.review-testing/interface-implementation/acousticHelmholtzFoam-baseline`.

## Geometric area fractions

Each valid adjacent PLIC plane supplies its own cut-face area. Two available
fractions are averaged; one is used directly. A valid zero or full fraction
still counts. The calculation is independent of `phi`.

Faces are triangulated using a canonical vertex ordering and the first valid
vertex fan. This fixes the triangulation under processor-face reversal.
Clipping uses the existing `cutFaceAdvect` static geometry operation. For
warped faces, scalar coverage uses triangle surface areas; the flux retains
oriented phase-area vectors whose sum is the original face-area vector.
Faces without a valid nonoverlapping vertex fan are rejected.

Positive signed distance from the PLIC plane is the liquid phase. Missing
planes in mixed cells are errors; pure-cell classification uses the selected
reconstruction scheme's `surfCellTol`. Opposite pure cells generate a
coincident-interface record. Their scalar fraction retains the original
interpolation fallback for the geometry-only comparison.

Processor neighbours exchange both plane descriptions before cutting.
Each copy computes the same unsigned average, and the code checks their
agreement to `1e-12`. Counts printed in logs include processor faces on both
ranks.

## Experimental transmission discretization

This opt-in prototype is retained for research and is excluded from the submission
configuration. To reproduce it, use `examples/transmission`; the numerical controls
are `maxStencilRings 3`, `svdRelativeTolerance 1e-12`, and `maxConditionNumber 1e10`.


For each candidate plane, pressure is reconstructed as

\[
P_j(x)=p_\Gamma+g_t\cdot(x-x_\Gamma)
       +\rho_jq_\Gamma n\cdot(x-x_\Gamma),\qquad g_t\cdot n=0.
\]

Pressure and inverse-density normal derivatives are continuous at that plane.
Owner and neighbour pressures are exact constraints; additional cell-centre
values determine the remaining coefficients by weighted least squares.
Coordinates and density are scaled before the nullspace/SVD calculation.
The stencil expands only as needed for the fit, up to the configured limit.

Each candidate's complete liquid-plus-gas face flux is calculated, then the
candidate fluxes are averaged. Pressure values and gradients are not averaged
independently of the geometric cuts. Elimination yields weights on existing
cell pressures; the solver keeps one complex pressure per cell.

A processor face has one builder. Its operator is shared with opposite
orientation on the neighbouring rank. Assembly conserves the inverse-density
pressure-gradient flux before applying the existing density row scaling.
The same real coefficients enter both pressure-component diagonal blocks.

The original physical Laplacian coefficient is zeroed on replaced faces.
This removes its implicit term and explicit nonorthogonal source there.
PML coupling, reaction quadrature, density mixing and existing velocity/force
postprocessing are unchanged. Geometry and reconstruction weights are cached
for the fixed mesh.

### Supported geometry and limitations

The implementation supports internal and processor faces, Cartesian 1D/2D/3D
meshes, and axisymmetric wedge meshes. Empty and symmetry faces retain their
boundary treatment. Mixed-cell support must stay away from nonsymmetry
physical boundaries and nonzero PML coefficients. Moving meshes and new-mode
cyclic/AMI/overset configurations are rejected.

A rank-deficient or excessively ill-conditioned fit is an error. There is no
automatic substitution of the legacy method.

**The transmission option is experimental.** Exact affine transmission and
conservative assembly do not guarantee accurate finite-frequency solutions.
In particular, the implemented method has a large error for a quarter-cell
shifted, high-contrast layered test at 20 kHz. An independent assembly
reproduces that error. The retained reaction quadrature and existing velocity
postprocessing must be considered when interpreting results. See
`testsuite/acousticInterface/VERIFICATION.md` for measured results and failed
acceptance gates.

## Diagnostics and verification

With `writeDiagnostics true`, `alphaf` and per-rank TSV files are written under
`postProcessing/acousticInterface/<time>/`:

- `geometry.tsv`: candidate fractions, oriented phase areas and PLIC planes.
- `operators.tsv`: global cell indices, flux weights, stencil rings and conditioning.
- `flux.tsv`: integrated inverse-density pressure-gradient flux for each component.

These fluxes are defined by the integral of `grad(P)/rho`; they are not
radiation-force or mass-flow outputs. Existing `Ure`, `Uim`, `pr` and `momFlux`
retain their original definitions.

The solver checks PETSc return codes and the convergence reason, rejects
nonfinite results, and reports the true assembled residual and backward error.
It also reports pressure changes between nonorthogonal corrections. The test
drivers check the latter when multiple corrections are requested.

Run the compiled and integration checks with Python, NumPy, SciPy and pytest:

```sh
python3 -m pytest --confcutdir=testsuite/acousticInterface \
    testsuite/acousticInterface/test_invariants.py
```

Run controlled comparisons without rebuilding:

```sh
python3 testsuite/acousticInterface/compare.py \
    --output /absolute/new/results-directory --ranks 1,2,4

python3 testsuite/acousticInterface/advanced.py \
    --output /absolute/new/sphere-results --scenario sphere \
    --sizes 16,24,32,48 --ranks 1,4

python3 testsuite/acousticInterface/wedge.py \
    --output /absolute/new/wedge-results --ranks 1,4
```

The wedge driver additionally requires cfMesh's `cartesian2DMesh`.
Each output directory must be new. Prepared meshes and initial fields are
copied for each method, and method settings are applied after preparation.
The drivers retain source/input/binary hashes and raw cases, and reject
nonfinite or incomplete data. `advanced.py` and `wedge.py` retain failed runs
in their result tables and return a nonzero exit status if any run fails.
A finished study must not be mistaken for passed accuracy or stability gates.
Python dependencies are listed in `testsuite/acousticInterface/requirements.txt`.
The analytical sphere driver also supports `--scenario sphereWedge`.
For oblique meshes, use `--mesh-style regular` or `--mesh-style perturbed`.

Lichtenberg build, study and wedge submission scripts are supplied alongside
the tests. They use isolated binaries, the existing PETSc installation, and
scratch output. Record the scheduler job IDs when submitting them.

The cluster wedge script expects an actual cfMesh source snapshot in
`.review-testing/interface-implementation/cfmesh-source`, including
`meshLibrary/Make/files`. The v2512 installation on Lichtenberg has an empty
cfMesh placeholder; the verification used a copy of the local v2606 bundled
cfMesh sources and built them against v2512 in isolated output directories.

Additional tools audit the actual matrix (`audit_operator.py`), shared
processor weights (`audit_parallel.py`), and retained reflection/transmission
phases and velocity diagnostics (`supplement_layered.py`). The verification
report records the test scope, scientific failures and scheduler job IDs.
