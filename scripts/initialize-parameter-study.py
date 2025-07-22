#!/usr/bin/env python3

from argparse import ArgumentParser
from subprocess import run
import os
import sys
import time


app_description = """ Initialize all test cases that match a given pattern.
    Initialization consists of the following steps:
        1. Create a mesh using the prescribed meshing application.\n
        2. Initialize required fields.\n
        3. Optional: perform domain decomposition for parallel simulation.\n
"""

def find_case_directories(pattern):
    """
    Find all directories whose name starts with 'pattern'.
    """
    case_dirs = [case_dir for case_dir in os.listdir(os.curdir) if os.path.isdir(case_dir) and \
                  case_dir.startswith(pattern)] 

    case_dirs.sort()

    return case_dirs

def submit_slurm_job(command, log_name=""):
    # partition = "test30m"
    account = "special00005"
    ntasks = "1"
    mem_per_cpu = "180G" # "18000"
    time_limit = "00:30:00" # Read: hh:mm:ss, h-> hours, m-> minutes, s-> seconds
    slurm_command = ("srun --account %s --ntasks %s --mem-per-cpu %s "
                     "--time %s "
                     "--job-name=%s --output=log.%s %s >/dev/null 2>&1 &" % 
                     (account, ntasks, mem_per_cpu,
                      time_limit,
                      log_name, log_name, command))

    run(slurm_command, shell=True)

    # Force wait after submitting job to avoid overloading the SLURM manager
    time.sleep(1)

def execute_init_step(command, use_slurm, log_name=""):
    if use_slurm:
        submit_slurm_job(command, log_name=log_name)
    else:
        run(command, shell=True)


#---- Command line arguments ----------------------------------------------
parser = ArgumentParser(description=app_description)

parser.add_argument("directory_pattern",
                    help="Pattern matching those directory names which are to be initialized, e.g. myAwesomeStudy_000")
parser.add_argument("-m", "--meshing-application",
                    help="Use this application for mesh creation, e.g. blockMesh.\nDefault: blockMesh",
                    required=True,
                    dest="mesh_app")
parser.add_argument("-f", "--field-init-script",
                    help="Name of the script to use for field initialization.\nDefault: None",
                    default="",
                    dest="init_script")
parser.add_argument("-j", "--job-mode",
                    help="Submit meshing and preprocessing as SLURM jobs",
                    action="store_true",
                    default=False,
                    dest="job_mode")
parser.add_argument("-par", "--parallel",
                    help="Use domain decomposition to setup cases for parallel execution.",
                    action="store_true",
                    default=False,
                    dest="parallel")

#--------------------------------------------------------------------------

if __name__ == '__main__':

    args = parser.parse_args(sys.argv[1:])

    if not args.init_script:
        print("No init script specified via --field-init-script. Not initializing fields.")


    case_dirs = find_case_directories(args.directory_pattern)
    
    case_count = 0
    for case in case_dirs:
        pwd = os.getcwd()
        os.chdir(case)

        case_count = case_count + 1
        print("(%s/%s) Initializing %s ..." % (str(case_count), str(len(case_dirs)), case))

        # Meshing
        # execute_init_step(args.mesh_app, args.job_mode, log_name=args.mesh_app)

        # Field initialization
        if args.init_script:
            # Note: the template replacement process invoked in 'argo-create-parameter-study.py'
            #       removes the She-Bang from the init script. So we need to use the file extension
            #       to detect the suitable interpreter (TT)
            script_env = ""
            if args.init_script.endswith('.py'):
                script_env = "python3"
            elif args.init_script.endswith('.sh'):
                script_env = "bash"
            else:
                sys.exit("Error: unsupported script format", args.init_script.rsplit('.')[0], ". Use either .sh for BASH or .py for Python3.")
            execute_init_step(script_env + " " + args.init_script, args.job_mode, log_name=args.init_script)

        # Domain decomposition
        if args.parallel:
            execute_init_step("decomposePar -force", args.job_mode, log_name="decomposePar")

        os.chdir(pwd)
