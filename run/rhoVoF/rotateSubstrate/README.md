# DropImpactRotatingPlate

This collection of cases was generated as part of the master thesis " Computational Fluid Dynamics and Experimental Analysis of Droplet Impact Dynamics and Heat Transfer Mechanisms on Rotating Plates" by Luis Meister.

## Contents

Contained within this directory are the 4 stages of case development namely:
	
* **reducedCase**: Initial case after parameterization script setup was lifted.
* **reducedCaseTestPLICaRDF**: Initial case adapted with function objects for PLIC and visRDF outputs, and parallelization setup for high cell count runs.
* **reducedCaseCHTtest**: Addition of Conjugate Heat Transfer to the adapted case, with imperfect modelling of the substrate rotation.
* **reducedCaseCHTtestRealVeloSR**: Exploratory setup with higher impact velocity and new modelling approach for substrate rotation.

Furthermore a selection of ParaView pipelines to analyse the results of the simulation:
	
* **ParaViewPipeline**

and a sbatch script file to replicate the start up process on the Lichtenberg supercomputer

* **SBATCHscriptForCluster**

## Requirements

To run this setup OpenFOAM v2412 is required, along with the compatible ThirdParty directory

The ParaView version for the pipeline is ParaView-5.13.2-MPI-Linux-Python3.10-x86_64

Additionaly the TwoPhaseFlow version set up in this branch should be used, as it was adapted to function with OF v2412.

## Installation

To build OpenFOAM v2412 with it's ThirdParty directory follow these steps

```bash
	cd $HOME && mkdir OpenFOAM && cd OpenFOAM	#creates an OpenFOAM folder to work in
	git clone https://develop.openfoam.com/Development/ThirdParty-common.git	#downloads/mirrors the third party dict
	mv ThirdParty-common ThirdParty-v2412		#renames the download
	cd ThirdParty-v2412 && git checkout v2412        #goes into dict and sets the thirdparty repo from latest version to v2412

	cd $HOME/OpenFOAM
	git clone https://develop.openfoam.com/Development/openfoam.git			#downloads/mirrors the OpenFOAM dict
	mv openfoam OpenFOAM-v2412			#renames the download
	cd OpenFOAM-v2412 && git checkout tag OpenFOAM-v2412				# goes into dict and sets the OpenFOAM repo from latest version to v2412
	source etc/bashrc				#sources environment data for OpenFOAM, tells pc where files are

	cd $HOME/OpenFOAM/ThirdParty-v2412		#goes back to Thirdparty
	./Allwmake					#should build thirdparty

	cd $HOME/OpenFOAM/OpenFOAM-v2412		#goes back to main OpenFOAM
	./Allwmake					#should build OpenFOAM
```



To install further needed ThirdParty libraries like kahip or scotch, clone their respective online repository into your local ThirdParty folder, rename it according to the folder name listed by the ./Allwmake process during building and rebuild ThirdParty, afterwards rebuild OpenFOAM. DO NOT SKIP THIS STEP, YOU'LL NEED THESE MODULES FOR PARALLEL RUNS!



To install TwoPhaseFlow run the following

```bash

	cd $HOME					#go to your main directory
	git clone https://github.com/tmaric/TwoPhaseFlow.git				#download/mirror TwoPhaseFlow
	cd TwoPhaseFlow && git checkout feature/density-ratio-DropImpactRotatingPlate	#go into folder and go to this branch
	cd modules/multiDimAMR && ./Allwmake						#go into AMR folder and build
	cd ../.. && ./Allwmake								#build TPF with this setup
	
	cd $HOME/OpenFOAM/OpenFOAM-v2412 && ./Allwmake	#rebuild OpenFOAM
```

## Run a case

To run a case go into a case directory, ensure that it is setup for the number of cores at your disposal (e.g. check system/balanceParDict and system/decomposeParDict) and execute the next steps in the following order:

```bash 
	source $HOME/OpenFOAM/OpenFOAM-v2412/etc/bashrc		#load OpenFOAM env
	source $HOME/TwoPhaseFlow/scripts/bashrc		#load env for TwoPhaseFlow PLIC and visRDF
	./Allclean			#Cleans up the case from any leftovers of previous simulation runs, use with care
	./Allrun.pre OR ./AllrunCHT	#execute all pre steps to set the case up for simulation
	./Allrun			#Only contains the command to run the case with the selected solver and with the set number of cores
```
ESPECIALLY FOR THE LATER CASES THE ALLRUN FILE IS OBSOLETE AND THE SIMULATION START IS PERFORMED DIRECTLY FROM SBATCH FILE!

A sbatch file is only used to execute the cases on a cluster, if run locally just execute it from an Allrun file (which might need to be adapted to the specific use case)

Finally in the provided sbatch file are also environmental source calls for OpenFOAM and TwoPhaseFlow which NEED to be included into the local workflow to run the cases properly (an example of this can be found in the previous code sample).

The source calls need to be executed every time you restart your terminal!

Attention for the cases should be placed on:

* **0.orig/U**	In here the RPM for the case can be set.
* **system/setAlphaFieldDict** Depending on the set RPM the placement of the droplet should be adapted, since the experimental data varies in it's resolution.

Further details on the other files can be either found in the master thesis or directly in the comment files.

Sometimes the compilation of the included Function Objects for the case **reducedCaseTestPLICaRDF** fails while starting the run for a simulation. Try deleting the "PolyMesh" file in the "constant" directory and restart the case, the folder should be automatically replaced with a new one and the Function Object compiled properly.

## Parallelization

Of the provided setups all apart from **reducedCase** are capable of parallel running the case on multiple CPU cores. You need to ensure that the value for "fileHandler" in controlDict is "uncollated" to generate the correct file setup for the respective amount of CPUs. Furthermore, the files "decompoaseParDict" and "balanceParDict" need to be adapted to the respective amount of CPU cores available.

All commands aimed at running the case in parallel need to be executed with "mpirun -np 8" before them (where 8 can be replace by any arbitrary number of available cores).
For example you would use "mpirun -np 1 ./Allrun.pre" to set up a non-CHT case for a parallel run and afterward execute it with "mpirun -np 96 $(getApplication) -parallel > logfile 2>&1". 
Further examples can be found in the file "SBATCHscriptForCluster".

## ParaView workflow

In ParaView the simulation results can be analysed by opening a ParaView-state file located within "ParaViewPipeline". Here the case file is currently selected to create the result data for the thesis but any other case can be loaded by modifying the currently shown base case. 

It should be ensured that the ParaView versions match especially if working on a cluster. The ParaView binaries can be found here: https://www.paraview.org/files/

In each ParaView state file at any given moment a sceenshot can be taken with the program. To automate this process and to create for example a sequence of overlay images a state file can be saved not only in ParaView format .pvsm but also in python format .py. 

This python script file will recreate the same state saved and can be edited in any way. 
Several example scripts of this type with adaptations for automatic sequencing can be found in the ParaView folder as well. 

To create a new sequencing script a new state can be saved in .py format and edited with:

```bash
	# get animation scene
	animationScene1 = GetAnimationScene()
	animationScene1.UpdateAnimationUsingDataTimeSteps()

	tk = GetTimeKeeper()
	tsteps = tk.TimestepValues
	#print(tsteps)
	i=0
	p=180  #(luis) start frame needs to be matched manually
	print(tk)
	print(tsteps)
	for t in tsteps:
	    tk.Time=tsteps[i]   #(luis) changes the timestep
	    rpm800_Layer2_Set40180 = CreateTexture(filename='/FILEPATH_EXPERIMENTAL_PICS/BW800RPM/Rpm800_Layer2_Set40'+str(p)+'.png')
```
and

```bash
	print(i)
	SaveScreenshot("/FILEPATH_OUTPUT_FOLDER/HighCellvisRDFContour"+str(i)+".png")
	i=i+1
	p=p+5
```
at the respective places in the scripts, thus creating a python for-loop over a given sequence.

Finally the example or newly created scripts should be executed with:

```bash
	pvpython --force-offscreen-rendering FILENAME.py --mesa
```
which allows for rendering without a screen, which is especially useful for a cluster application.
