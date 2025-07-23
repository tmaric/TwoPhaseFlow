#!/bin/bash
#SBATCH -J spreading
# #SBATCH -o jetCrossFlow.out.%j
# #SBATCH -e jetCrossFlow.err.%j
#SBATCH -n 100
# #SBATCH -c 1
#SBATCH --mem-per-cpu=8000
#SBATCH -t 06-23:00:00
#SBATCH -A special00004
# #SBATCH --mail-user=jliu.fb16@gmail.com
# #SBATCH --mail-type=all

mpirun -np 100 interFlow   -nonOrthogControl PIMPLE -parallel
