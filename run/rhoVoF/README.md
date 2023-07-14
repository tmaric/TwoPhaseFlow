# test cases for cyclic BC debugging and interIsoRhoFoam

benchmark cases for verification

NOTE: the cyclic BC for interIsoFoam has not solved util OpenFOAM-V2212. To not pollute the origin codes in OpenFOAM, we modified the BC handling only in TwoPhaseFlow.
We configured the cases so that interFlow works same as interIsoFoam to avoid the cyclic BC problem of default interIsoFoam in OpenFOAM.

## cyclicBC_debug_case2D 

run case by

```bash
    ./Allrun
```

## jetCrossFlow

install cfMesh before running case.

preprocess case by 

```bash
    ./initialCase.sh
```

run case on cluster with SLURM

use interIsoFoam
```bash
    sbatch isoAdv_slurm
```

use interIsoRhoFoam
```bash
    sbatch isoRho_slurm
```

run case on local workstation

use interIsoFoam 
```bash
    foamJob -log-app -parallel interFlow 
```

use interIsoRhoFoam
```bash
   foamJob -log-app -parallel interIsoRhoFoam -tScheme Euler 
```

### mixingLayer2D, risingBubble3D, translatingDroplet3D, translatingDropletInQuiescentFluid3D

preprocess case by

```bash
    ./run_case_serial.sh [PARAMETER_NAME].parameter [PREFIX]
```

run case on cluster with SLURM

use interIsoFoam
```bash
    bulkval [CASES_NAME_COMMON_PATTERN] "sbatch ../isoAdv.sbatch"
```

use interIsoRhoFoam
```bash
    bulkval [CASES_NAME_COMMON_PATTERN] "sbatch ../isoRho.sbatch"
```

run case on local workstation

use interIsoFoam  
```bash
    bulkval [CASES_NAME_COMMON_PATTERN] "foamJob -log-app interFlow" 
```

use interIsoRhoFoam
```bash
   bulkval [CASES_NAME_COMMON_PATTERN]  "foamJob -log-app interIsoRhoFoam -tScheme Euler"
```

postprocess

run .ipynb file in each case. NOTE: adapt the name of first case in dataAgglomeration

