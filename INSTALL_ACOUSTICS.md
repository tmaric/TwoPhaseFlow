# Acoustic Solver Installation

This guide sets up the frequency-domain, time-domain, and coupled acoustic
solvers from a fresh clone. Commands assume a Linux workstation and Bash.

## Supported Environment

- OpenFOAM-v2606 from openfoam.com is the canonical target.
- OpenFOAM.org releases are not supported.
- The dependency revisions and acoustic smoke cases below have been reproduced
  with OpenFOAM-v2606.
- PETSc is built locally inside the repository with MPI, MUMPS, Hypre,
  SuperLU_DIST, ScaLAPACK, METIS, and ParMETIS.

## 1. Install And Source OpenFOAM

Install OpenFOAM-v2606 using the instructions supplied by openfoam.com. Source
the installation before configuring this repository:

```bash
source /path/to/OpenFOAM-v2606/etc/bashrc
echo "$WM_PROJECT_VERSION"
```

The printed version must be `v2606`.

## 2. Clone The Acoustic Branch

```bash
git clone --branch feature/acousticLevitation \
    https://github.com/tmaric/TwoPhaseFlow.git
cd TwoPhaseFlow
source ./scripts/bashrc
```

In every new terminal, source OpenFOAM first and this repository second:

```bash
source /path/to/OpenFOAM-v2606/etc/bashrc
source /path/to/TwoPhaseFlow/scripts/bashrc
```

## 3. Build PETSc And petsc4Foam

The dependency build can install its Ubuntu build prerequisites and may take
considerable time:

```bash
./petsc_build.sh \
    --install-deps \
    --openfoam-bashrc "$WM_PROJECT_DIR/etc/bashrc" \
    -j "$(nproc)"
```

The default dependency revisions are pinned for repeatability:

```text
PETSc:      0eed7b4d3dc0b28ce2d9ef5959622812f5883345
petsc4Foam: d72de946b0d893a7bb9b8e052cec2dad5bc1af73
```

Alternative commits or tags can be selected with `--petsc-ref` and
`--petsc4foam-ref`. Local dependency trees are created as `petsc/` and
`external/petsc4Foam/`; do not commit them.

Reload the project environment after the dependency build:

```bash
source ./scripts/bashrc
```

## 4. Build TwoPhaseFlow

The top-level build compiles the core libraries, `setPMLFields`, and the three
acoustic solvers:

```bash
./Allwmake
./apps/utilities/preProcessing/setAlphaField/Allwmake
```

The second command builds the repository's geometric `setAlphaField` utility,
which is required by the two-phase acoustic cases.

Build the cfMesh plugin supplied with OpenFOAM when using cases based on
`cartesian2DMesh`:

```bash
"$WM_PROJECT_DIR/plugins/cfmesh/Allwmake" -j "$(nproc)"
```

## 5. Check The Environment

Run the non-destructive checker after compilation:

```bash
./scripts/check-acoustic-environment.sh
```

It checks the OpenFOAM version, repository environment, PETSc installation,
required commands, solver executables, and unresolved shared libraries.

## 6. Python And Optional Tools

The complete verification-case workflows use NumPy, SciPy, Matplotlib, and
VTK for postprocessing:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy scipy matplotlib vtk oftest
```

Some legacy meshes require the `gmsh` command. Install it through the operating
system or the official Gmsh distribution. ParaView is optional but recommended
for inspecting meshes and fields.

Initialize the AMR module only when AMR is required:

```bash
git submodule update --init --recursive
./modules/multiDimAMR/Allwmake
```

## 7. Smoke Tests

Inspect each case's `Allrun`, `caseParams.sh`, `controlDict`, and
`decomposeParDict` before execution.

Frequency-domain smoke test:

```bash
cd run/acousticTests/FrequencyDomainTests/layeredInterface1D
NX=400 ./Allrun
tail -30 log.acousticHelmholtzFoam
```

Time-domain smoke test:

```bash
cd run/acousticTests/TimeDomainTests/layeredInterface1D
NX=400 ./Allrun
tail -30 log.acousticWaveFoam
```

A smoke test passes when the build and solver finish, the log reaches `End`,
and the expected pressure fields are finite. It does not establish analytical
accuracy, convergence, or publication evidence.

Coupled `interFALFlow` smoke test:

```bash
cd run/acousticTests/FrequencyDomainTests/Andrade2019/WedgeLevitatorDynamicDrop
CFMESH_MAX_CELL_SIZE=1e-3 \
END_TIME=1e-5 \
WRITE_INTERVAL=1e-5 \
./Allrun
tail -30 log.interFALFlow
```

The runtime overrides are applied only for the run. `Allrun` restores the
production `controlDict` when it exits.

## 8. Troubleshooting

Confirm the environment in the current shell:

```bash
echo "$WM_PROJECT_VERSION"
echo "$TPF_PROJECT_DIR"
echo "$PETSC_DIR/$PETSC_ARCH"
```

If a solver cannot find PETSc at runtime, source the environments again in the
required order:

```bash
source /path/to/OpenFOAM-v2606/etc/bashrc
source /path/to/TwoPhaseFlow/scripts/bashrc
```

Then rerun `./scripts/check-acoustic-environment.sh`.
