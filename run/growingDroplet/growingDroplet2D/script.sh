#!/bin/bash -l

#allocate 1 node for 24 hour on LB2 partition
#SBATCH -A special00005 
#SBATCH -n 144 
#SBATCH --mem-per-cpu=3000 
#SBATCH -t 23:30:00     
#SBATCH --error=error.err
#SBATCH -J zeroG 

module load git
module load gcc/9.2.0
module load openmpi/4.0.3
module load boost
module load eigen/3.3.7
module load cmake
module load python/3.7.4
module load gnuplot
module load texlive
source /work/projects/special00005/B01/OpenFOAM-v2112/etc/bashrc

srun -n 144 interFlow -parallel -fileHandler collated > log
srun -n 144 foamToVTK -parallel  > log.vtk

#reconstructPar -fileHandler collated
