# How to run jetCrossFlow case with AMR

## Build multiDimAMR module in TwoPhaseFlow 

```bash
cd TwoPhaseFlow
./build-multiDIMAMR.sh
```
The module is ported to the OpenFOAM+ version v1812 and v2006,v2012. For new version, code adapation is necessary.

## Install/build kahip

```bash
git clone https://develop.openfoam.com/Development/ThirdParty-common.git
mv ThirdParty-common ThirdParty
mv ThirdParty $WM_PROJECT_DIR
cd $WM_PROJECT_DIR/ThirdParty
./Allwmake
cd $FOAM_SRC/parallel/decompose
./Allwmake
```

##Run simulation
```bash
cd jetCrossFlow_AMRkahip/
./initialCase.sh
mpirun -np 100 interFlow -parallel 2>&1 | tee log.interFlow
```



