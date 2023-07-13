#!/bin/bash
# #SBATCH -p test7d
#SBATCH -A special00004
#SBATCH -J interIsoFoam
#SBATCH -n 1
#SBATCH --mem-per-cpu=8000
#SBATCH -t 04-23:00:00

srun interFlow #interIsoFoam
