import sys
import shlex
import time
import os
import argparse
import yaml

parser=argparse.ArgumentParser(epilog=('#### Any other arguments passed will be passed as they are to prosculpt. If including output_dir, please include it first ####\n'))

parser.add_argument('num_tasks', help='Number of concurrent jobs to run.')
parser.add_argument('name', help='Desired task name')
parser.add_argument('-d', '--dry-run', action="store_true", help="Print command but do not run.")
args=parser.parse_known_args() #This allows us to pass any other arguments further to prosculpt


WIDTH = 80
SYMBOL = "⁝"
print("".center(WIDTH, SYMBOL))
print(" prosculpt2slurm ".center(WIDTH, SYMBOL))
print("".center(WIDTH, SYMBOL))
print("Please make sure your first argument was the number of concurrent jobs you want to run and the second argument was the desired name of the slurm task. After that, add the arguments that will be passed as-is to the _merged script.")
print("".center(WIDTH, SYMBOL))

slurm_runner_path= os.path.dirname(os.path.realpath(__file__))
current_path= os.getcwd()

n = int(args[0].num_tasks)
task_name = args[0].name

if task_name[:11] == "output_dir=":
    print(f"!!! Are you sure you want to name your slurm job {task_name}? Anyway, continuing ...")

out_command_file = f"ps2slurm_{task_name}_{int(time.time())}.txt"

with open(out_command_file, 'w') as f:
    for i in range(1,n+1):
        let_argsv = args[1].copy()

        output_argn=-1
        for argn, arg in enumerate(let_argsv):
            if "output_dir" in arg: output_argn=argn
            if "-cd" in arg: config_file_folder_argn=argn+1
            if "-cn" in arg: config_file_name_argn=argn+1

        if len(let_argsv)>0:
            if output_argn != -1: #If output_dir is in the let_argsv
                if let_argsv[output_argn][-1:]!="/":
                    let_argsv[output_argn]+="/"
                let_argsv[output_argn] += f"{i:02d}"
            else:  #if not, it is in the yaml file. Read it and modify it if there are more than 1 task
                #with open(os.path.join(slurm_runner_path,let_argsv[config_file_folder_argn],(let_argsv[config_file_name_argn])+".yaml"), "r") as yaml_file:
                with open(os.path.join(current_path,let_argsv[config_file_folder_argn],(let_argsv[config_file_name_argn])+".yaml"), "r") as yaml_file:
                    data = yaml.safe_load(yaml_file)
                output_dir=data["output_dir"]
                if output_dir[-1:]!="/":
                    output_dir+="/"
                output_dir += f"{i:02d}"
                let_argsv.append("++output_dir="+output_dir) #This will override the output present in the yaml file
                #print("+output_dir="+output_dir)

        cmdline = " ".join(map(shlex.quote, let_argsv)) #join all arguments passed that aren't number of tasks or task name
        line = f"""python {slurm_runner_path}/rfdiff_mpnn_af2_merged.py {cmdline}"""
        print(line, file=f)

print(f"Slurm command can be found in {out_command_file}")

if not args[0].dry_run:
    #Exclude is there because that node was working incredibly slow
    os.system(f"export GROUP_SIZE=1; sbatch --exclude=compute-0-11 --partition=gpu --gres=gpu:A40:1 --ntasks=1 --cpus-per-task=2 --output slurm-%A_%a_{task_name}.out --error slurm-%A_%a_{task_name}.out -J {task_name}  -a 1-{n} {slurm_runner_path}/wrapper_slurm_array_job_group.sh {out_command_file}")
    print(f"Job {task_name} has been submitted to slurm".center(WIDTH, SYMBOL))
else:
    print("Command wasn't run because --dry-run was active.")