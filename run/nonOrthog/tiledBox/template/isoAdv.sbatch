#!/bin/bash
#SBATCH -J mixingLayer
# #SBATCH -o jetCrossFlow.out.%j
# #SBATCH -e jetCrossFlow.err.%j
#SBATCH -n 1
# #SBATCH -c 1
#SBATCH --mem-per-cpu=8000
#SBATCH -t 23:00:00
#SBATCH -A special00004
# #SBATCH --mail-user=jliu.fb16@gmail.com
# #SBATCH --mail-type=all

interFlow #interIsoFoam
