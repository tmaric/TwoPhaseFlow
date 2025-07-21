# TwoPhaseFlow

The TwoPhaseFlow library adds new surface tension and phase change models to OpenFOAM and provides benchmark cases for verification.

## Documentation

The available models and solvers are documentated in the paper:

Scheufler, H., & Roenby, J. (2023). "TwoPhaseFlow: A Framework for Developing Two Phase Flow Solvers in OpenFOAM". OpenFOAM® Journal, 3, 200–224. https://doi.org/10.51560/ofj.v3.80

Extensions for handling high density ratios are documented in: 

[Liu, J., Tolle, T., Zuzio, D., Estivalèzes, J. L., Damian, S. M., & Marić, T. (2024). Inconsistencies in unstructured geometric volume-of-fluid methods for two-phase flows with high density ratios. Computers & Fluids, 281, 106375.](https://doi.org/10.1016/j.compfluid.2024.106375)

## Getting Started


### Prerequisites

Requires OpenFOAM-v1812 or later:

```
https://www.openfoam.com/download/release-history.php
```

Please checkout the appropriate branch to compile with later OpenFOAM version.  

OpenFOAM.org versions are not supported.

### Installing

```bash
    git clone https://github.com/DLR-RY/TwoPhaseFlow
    cd TwoPhaseFlow
    # To compile e.g. with OpenFOAM-v2506 checkout the appropriate branch with:
    # git checkout of2506
    ./Allwmake
    ./get-gmsh.sh # will install gmsh version 493 as gmshv493
    # for AMR
    git submodule update --init --recursive
    cd modules/multiDimAMR/
    ./Allwmake
```
### Running testsuite

Make sure that the desired OpenFOAM installation is sourced e.g. v2506 and that 
python is installed with a version >= 3.6 (miniconda is a great option, but anaconda works as well)

```bash
    python -m venv env # creats virtual python enviroments (optional step)
    pip install oftest

    py.test # runs the tests
    py.test --writeNSteps=1 run/ # test all testcases in run
```

## Authors

* **Henning Scheufler**

## Contributors

* **Tomislav Maric**
* **Tobias Tolle**
* **Jun Liu**

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

