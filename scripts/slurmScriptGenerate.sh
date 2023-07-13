#!/bin/bash

#--------------------------------------------------------------------
# Slurm script generator for running OpenFOAM cases on HPC clusters.
#
# Note:
#   Please put this file in the root directory of the case.

# Author:
#   Dezhi Dai (dezhi.dai@mavs.uta.edu), MAE Department, UTA.
#--------------------------------------------------------------------

#-Modify this part if needed-----------------------------------------
# Set the email, replace this by your own email address
email="jliu.fb16@gmail.com"
# Set the flow solver
flow_solver="interIsoRhoFoam"
# Set the No. of nodes requested and total mpi tasks
n_proces="1"
n_cores_perProc="32"
# Set the project to charge
proj="special00005"
# Set request time
time_request="12:00:00"
# Main memory in MByte for each cpu core
mem_cpu="1500"
#--------------------------------------------------------------------

# Obtain current case directory and rootcase name
case_dir=`pwd` 
case_name=`basename "$case_dir"`
# Set the file name of the SLURM job script
s_name="auto_slurm_job"

# Write options and commands to the Slurm script
# Header
echo "#!/bin/bash" > $s_name

# Slurm options 
# By default both standard output and standard error are directed 
# to a file of the name "slurm-%j.out", 
# where the "%j" is replaced with the job allocation number.
echo "#SBATCH -J $case_name" >> $s_name
echo "#SBATCH -o $case_name.out.%j" >> $s_name
echo "#SBATCH -e $case_name.err.%j" >> $s_name
echo "#SBATCH -n $n_proces" >> $s_name
echo "#SBATCH -c $n_cores_perProc" >> $s_name
echo "#SBATCH --mem-per-cpu=$mem_cpu" >> $s_name
echo "#SBATCH -t $time_request" >> $s_name
echo "#SBATCH -A $proj" >> $s_name
echo "#SBATCH --mail-user=$email" >> $s_name
echo "#SBATCH --mail-type=all" >> $s_name

# OpenFOAM commands JUN: need we parallelize meshing?
echo "decomposePar -case $case_dir -fileHandler collated" >> $s_name
echo "mpirun -np $n_cores_perProc $flow_solver -parallel -case $case_dir -fileHandler collated" >> $s_name
#echo "reconstructPar -case $case_dir -fileHandler collated" >> $s_name
echo "touch $case_name.foam" >> $s_name
