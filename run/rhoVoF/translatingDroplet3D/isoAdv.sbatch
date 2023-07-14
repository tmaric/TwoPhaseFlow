#!/bin/bash
# #SBATCH -p test24
#SBATCH -A special00005
#SBATCH -J interIsoFoam
#SBATCH -n 1
#SBATCH --mem-per-cpu=8000
#SBATCH -t 03-15:00:00

srun interFlow #interIsoFoam
