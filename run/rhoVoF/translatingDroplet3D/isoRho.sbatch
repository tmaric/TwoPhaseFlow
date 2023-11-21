#!/bin/bash
# #SBATCH -p test24
#SBATCH -A special00005
#SBATCH -J interIsoRhoFoam
#SBATCH -n 1
#SBATCH --mem-per-cpu=8000
#SBATCH -t 03-12:00:00

srun interIsoRhoFoam -tScheme Euler
