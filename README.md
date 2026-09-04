# TwoPhaseFlow

The TwoPhaseFlow library adds new surface tension and phase change models to OpenFOAM and provides benchmark cases for verification.

## Documentation

The available models and solvers are documentated in the paper:

Scheufler, H., & Roenby, J. (2023). "TwoPhaseFlow: A Framework for Developing Two Phase Flow Solvers in OpenFOAM". OpenFOAM® Journal, 3, 200–224. https://doi.org/10.51560/ofj.v3.80

## Getting Started


### Prerequisites

The acoustic-development branch uses OpenFOAM-v2606 from openfoam.com. Other
OpenFOAM distributions are not supported for this workflow:

```
https://www.openfoam.com/download/release-history.php
```

Use the branch matching the intended OpenFOAM release. OpenFOAM.org versions
are not supported.

### Installing

```bash
    git clone --branch feature/acousticLevitation \
        https://github.com/tmaric/TwoPhaseFlow.git
    cd TwoPhaseFlow
    source /path/to/OpenFOAM-v2606/etc/bashrc
    source ./scripts/bashrc
    ./petsc_build.sh --openfoam-bashrc "$WM_PROJECT_DIR/etc/bashrc"
    ./Allwmake
```

### Acoustic workflow

PETSc and MUMPS must be built before the top-level `./Allwmake`, since that
script now includes the PETSc-backed acoustic solvers. See
[INSTALL_ACOUSTICS.md](INSTALL_ACOUSTICS.md) for the complete fresh-clone
procedure, dependency revisions, environment checks, and smoke tests.

### Running testsuite

Make sure that OpenFOAM-v2606 is sourced and that Python 3 is installed.

```bash
    python3 -m venv env
    source env/bin/activate
    pip install oftest

    pytest
    pytest --writeNSteps=1 run/
```

## Authors

* **Henning Scheufler**

## Contributors

* **Chuanchao Xu**
* **Jun Liu**
* **Tomislav Maric**


### adaptive mesh refinement with multiple regions

AMR with multiple regions does not work in version of1812 but it is fixed in newer versions.


To fix this, apply the patch (assumes OpenFOAM is already source):

```bash
    cp  patches/multiRegionAMR.patch $WM_PROJECT_DIR
    cp  patches/tableBase.patch $WM_PROJECT_DIR
    cp  patches/surfaceFieldValue.patch $WM_PROJECT_DIR
    cd $WM_PROJECT_DIR
    git apply multiRegionAMR.patch
    git apply tableBase.patch
    git apply surfaceFieldValue.patch

```
details see:

https://develop.openfoam.com/Development/openfoam/-/issues/1676

https://develop.openfoam.com/Development/openfoam/-/issues/1753

## License

This project is licensed under the GPL v3 License - see the [LICENSE.md](LICENSE.md) file for details


## Running benchmarks

```bash
    ./get-gmsh.sh # install gmsh
    pip install casefoam

```

The run/benchmark cases are run with


```bash
    cd run/benchmark/phaseChange/suckingInterface/
    python genCases.py # generates the study based and template case (here StefanProblem)
    ./Allrun # runs all the created testcases
    python getResults.py # to see results
```

Alternatively, the runAll.sh can be executed in the folder.

Note:

Some cases use the slurm queuing system and call `sbatch Allrun_Slurm` in the Allrun script, so you might need to modify it in the template case.
