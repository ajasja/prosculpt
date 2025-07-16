import sys
import shlex
import time
import os
import argparse
import yaml

parser=argparse.ArgumentParser(epilog=('#### Any other arguments passed will be passed as they are to prosculpt. If including output_dir, please include it first ####\n**** Important: hydra config overrides should be put before -cd and -cn ****\nExample: python slurm_runner.py 1 multipassinpaintseq_throw  output_dir="Examples/Examples_out/multipass_inpaintseq_throw" +throw=1 +cycle=0 -cd Examples -cn multipass_inpaintseq'))

parser.add_argument('yaml_file', help='Prosculpt input file.')
parser.add_argument('-d', '--dry-run', action="store_true", help="Print command but do not run.")
args=parser.parse_known_args() #This allows us to pass any other arguments further to prosculpt

slurm_runner_path= os.path.dirname(os.path.realpath(__file__))
current_path= os.getcwd()


yaml_file_path = args[0].yaml_file
extra_args = args[1].copy()
yaml_file_dir, yaml_file_name = os.path.split(yaml_file_path)
yaml_file_name_no_extension = yaml_file_name.rsplit( ".", 1 )[ 0 ]

with open(yaml_file_path) as yaml_file:
    yaml_data = yaml.safe_load(yaml_file)

with open(os.path.join(slurm_runner_path,"config", "installation.yaml")) as installation_yaml_file:
    installation_yaml_data = yaml.safe_load(installation_yaml_file)

installation_slurm_yaml_data=installation_yaml_data["slurm"]
slurm_data=yaml_data["slurm"]
n=yaml_data["num_tasks"]

if "task_name" not in yaml_data:
    task_name="prosculpt_task"
else:
    task_name=yaml_data["task_name"]

out_command_file = f"ps2slurm_{task_name}_{int(time.time())}.txt"



with open(out_command_file, 'w') as f:
    for i in range(1,n+1):
        arguments=[]

        output_dir=""
        output_dir_in_args=False
        for argn, arg in enumerate(extra_args):
            if "output_dir" in arg:
                output_dir_in_args=True
                output_dir=arg
                if output_dir[-1:]!="/": #Add the task number
                    output_dir+="/" 
                output_dir += f"{i:02d}"
                arguments.append(output_dir) #This will override the output present in the yaml file with the one with the tasknumber

        if not output_dir_in_args:
            output_dir=yaml_data["output_dir"]
            if output_dir[-1:]!="/": #Add the task number
                output_dir+="/" 
            output_dir += f"{i:02d}"
            arguments.append("++output_dir="+output_dir) #This will override the output present in the yaml file with the one with the tasknumber

        arguments.append("-cd "+ yaml_file_dir) 
        arguments.append("-cn "+ yaml_file_name_no_extension) 
        
        #print("+output_dir="+output_dir)

        cmdline = " ".join(arguments) #join all arguments passed that aren't number of tasks or task name
        line = f"""python {slurm_runner_path}/rfdiff_mpnn_af2_merged.py {cmdline}"""
        print(line, file=f)

print(f"Slurm command can be found in {out_command_file}")

options_string=""
job_specific_keys=[]
if slurm_data is not None:
    for key, value in slurm_data.items(): 
        job_specific_keys.append(key)
        if key=="slurm_options_string":
            options_string+= f" {value}"
        else:
            if len(key)==1:
                options_string+= " -"
            else:
                options_string+= " --"
            options_string+= f"{key} {value}"

#now take the default values from installation.yaml if not present in job yaml
for key, value in installation_slurm_yaml_data.items(): 
    if key not in job_specific_keys:
        if len(key)==1:
            options_string+= " -"
        else:
            options_string+= " --"
        options_string+= f"{key} {value}"
        

if not args[0].dry_run:
    full_command= f"export GROUP_SIZE=1; sbatch -J {task_name} -a 1-{n} {options_string} {slurm_runner_path}/wrapper_slurm_array_job_group.sh {out_command_file}"
    print(f"Full command is: {full_command}")
    exit_code = os.system(full_command)
    if exit_code == 0:
        print(f"Job {task_name} has been submitted to slurm with code {exit_code}")
    else:
        print(f"Job submission failed with code {exit_code}")
else:
    print("Command wasn't run because --dry-run was active.")