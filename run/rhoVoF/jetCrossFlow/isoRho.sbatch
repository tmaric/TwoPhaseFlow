#!/bin/bash
#SBATCH -J jetCrossFlow
#SBATCH -n 32
#SBATCH -c 1
#SBATCH --mem-per-cpu=40000
#SBATCH -t 03-23:00:00
#SBATCH -A special00004

mpirun -np 32 interIsoRhoFoam -tScheme Euler -parallel #-fileHandler collated -parallel 
