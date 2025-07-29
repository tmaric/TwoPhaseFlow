#!/usr/bin/env python3

from argparse import ArgumentParser
import datetime 
import os
import shlex        # Parse shell arguments for solver, see
                    # https://docs.python.org/3/library/subprocess.html#subprocess.Popen
import subprocess
import sys
import time

import testReportCore as trc

app_description="""Execute a parameter study with a user-specified solver.
Can execute n variants in parallel.
"""

# The solution to control the number parallel running processes is taken / adapted from
# https://stackoverflow.com/questions/18123325/always-run-a-constant-number-of-subprocesses-in-parallel

nextVariant = 0
maxNumProcesses = 1
nVariants = 1
solver = ""
Processes = []

def start_new(variantList, n_mpi, solver_args):
    """ Start a new subprocess if there is work to do """
    global nextVariant
    global Processes

    args = shlex.split(solver_args)

    if nextVariant < nVariants:
        # Make the script wait for the last spawned subprocess.
        # This makes it more likely that the end time info
        # will be correct
        solver_command = [solver] + args + ["-case", variantList[nextVariant]]
        mpi_prefix = ["mpirun", "-n", str(n_mpi)]

        if nextVariant == nVariants-1:
            if n_mpi > 1:
                proc = subprocess.Popen(mpi_prefix + solver_command + ["-parallel"]).wait()
            else:
                proc = subprocess.Popen(solver_command).wait()
        else:
            if n_mpi > 1:
                proc = subprocess.Popen(mpi_prefix + solver_command + ["-parallel"])
            else:
                proc = subprocess.Popen(solver_command)

        print ("Started to process variant", variantList[nextVariant].rsplit('/',1)[1])
        nextVariant += 1
        Processes.append(proc)

def check_running(variantList, n_mpi, solver_args):
   """ Check any running processes and start new ones if there are spare slots."""
   global Processes
   global nextVariant

   for p in range(len(Processes)-1,0,-1): # Check the processes in reverse order
      if Processes[p].poll() is not None: # If the process hasn't finished will return None
         del Processes[p] # Remove from list - this is why we needed reverse order

   while (len(Processes) <= maxNumProcesses) and (nextVariant < nVariants): # More to do and some spare slots
      start_new(variantList, n_mpi, solver_args)


def main():
    
    global solver
    global maxNumProcesses
    global nVariants

    #---- Command line arguments ----------------------------------------------
    parser = ArgumentParser(description=app_description)
    parser.add_argument("solver",
                        help="Name of the solver to run the study with.")
    parser.add_argument("--solver-args",
                        help="Additional arguments to pass to the solver. This option is ignored if SLURM is used.\n Default: None",
                        type=str,
                        default="",
                        dest="solver_args")
    parser.add_argument("-d","--dir-pattern",
                        help="Pattern of the directories belonging to the study.",
                        type=str,
                        required=True,
                        dest="dir_pattern")
    parser.add_argument("-n","--num-processes",
                        help="The number of processes running in parallel.\n Default: 1",
                        type=int,
                        default=1,
                        dest="numprocesses")
    parser.add_argument("-j", "--job-mode",
                        help="Submit solver execution as SLURM jobs. Expects to find a SLURM script named solver_name.sbatch.",
                        action="store_true",
                        default=False,
                        dest="job_mode")
    parser.add_argument("-u","--use-mpi",
                        help="Call solver in parallel using the given number of processes.\n Has no effect if --job-mode is used.",
                        type=int,
                        default=0,
                        dest="num_mpi_procs")

    args = parser.parse_args()

    # Get list of variation folders
    studyDirectories = trc.pattern_cases(args.dir_pattern + "*")

    if not studyDirectories:
        print("Error: no directories found for study",args.dir_pattern,".")
        print("Exiting.")
        sys.exit()

    # Set the global parameters for parallel run
    solver = args.solver
    nVariants = len(studyDirectories)
    maxNumProcesses = args.numprocesses

    #--------------------------------------------------------------------------
    #   SLURM batch mode
    #--------------------------------------------------------------------------
    if args.job_mode:
        job_script = solver + ".sbatch";
        if not os.path.isfile(job_script):
            sys.exit("Error: no SLURM script named '" + job_script + "' found.\nExiting.")

        print("Submitting jobs to SLURM using job script", job_script)

        for index in range(len(studyDirectories)):
           pwd = os.getcwd()
           os.chdir(studyDirectories[index])

           print("\nSubmitting job", index+1, "/", len(studyDirectories))
           subprocess.run(["sbatch", "../" + job_script])
           time.sleep(1)

           os.chdir(pwd)

    else:
    #--------------------------------------------------------------------------
    #   Direct execution
    #--------------------------------------------------------------------------
        # User Info
        print("------------ Run Info ----------------")
        print("Running",args.dir_pattern,"using the solver",solver,"with",maxNumProcesses,
                "processes.")
        if args.num_mpi_procs > 1:
            print("Running with MPI using", args.num_mpi_procs, "mpi processes each.")
        print("Start_time: ",datetime.datetime.now().time())

        # Spawn processes
        check_running(studyDirectories, args.num_mpi_procs, args.solver_args) # This will start the max processes running
        while (nextVariant < nVariants): # Some thing still going on.
            time.sleep(0.1) # You may wish to change the time for this
            check_running(studyDirectories, args.num_mpi_procs, args.solver_args)

        print("End_time: ",datetime.datetime.now().time())
        print("Done!")

if __name__ == "__main__":
    main()
