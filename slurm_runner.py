import sys
import shlex
import time
import os
import argparse
import subprocess
import yaml

parser=argparse.ArgumentParser(epilog=('#### Any other arguments passed will be passed as they are to prosculpt. If including output_dir, please include it first ####\n**** Important: hydra config overrides should be put before -cd and -cn ****\nExample: python slurm_runner.py 1 multipassinpaintseq_throw  output_dir="Examples/Examples_out/multipass_inpaintseq_throw" +throw=1 +cycle=0 -cd Examples -cn multipass_inpaintseq'))

parser.add_argument('yaml_file', help='Prosculpt input file.')
parser.add_argument('-d', '--dry-run', action="store_true", help="Print command but do not run.")
parser.add_argument(
    '--allow-custom-log-path',
    action="store_true",
    help=(
        "By default, sbatch's own -o/-e (slurm's stdout/stderr log files) are forced to "
        "<output_dir>/logs/, regardless of any slurm.output/slurm.error set in the job yaml or in "
        "installation.yaml - so a job's slurm log always lives next to its own results, not wherever "
        "a job/site config happened to point it. Pass this to opt back into respecting slurm.output/"
        "slurm.error instead (the pre-existing behavior)."
    ),
)
args=parser.parse_known_args() #This allows us to pass any other arguments further to prosculpt

slurm_runner_path= os.path.dirname(os.path.realpath(__file__))
current_path= os.getcwd()


yaml_file_path = args[0].yaml_file
yaml_file_path=os.path.abspath(yaml_file_path)

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

out_command_file = f"ps2slurm_{task_name}_{int(time.time_ns())}.txt"

# Base output_dir - the same value the per-task loop below starts from
# before appending each task's own "01"/"02"/... suffix - used, unless
# --allow-custom-log-path is passed, as where sbatch's own -o/-e (slurm's
# stdout/stderr log files, as opposed to whatever prosculpt itself later
# logs to) get forced to instead. This can only be the *base* output_dir,
# not each task's own numbered one: all n tasks share a single sbatch
# array submission (one -o/-e template for the whole array), and SLURM's
# own %a substitution isn't zero-padded to match "01"/"02"/..., so there's
# no way for a shared template to land in each task's own numbered
# directory - the parent all of them share is the one thing that's
# actually common across the whole array.
base_output_dir = None
for arg in extra_args:
    if "output_dir" in arg:
        base_output_dir = arg.split("=", 1)[-1]
if base_output_dir is None:
    base_output_dir = yaml_data["output_dir"]

forced_log_dir = None
if not args[0].allow_custom_log_path:
    forced_log_dir = os.path.join(base_output_dir, "logs")
    # sbatch does not create the directory portion of -o/-e itself - if it
    # doesn't already exist when the job is submitted, the job fails
    # immediately with no log written anywhere at all (see
    # job_staging.ensure_logs_dir() in the dashboard for the same problem
    # solved the same way, for the cwd-relative logs/ this replaces).
    os.makedirs(forced_log_dir, exist_ok=True)


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

        arguments.append(f"-cd '{yaml_file_dir}'") 
        arguments.append(f"-cn '{yaml_file_name_no_extension}'") 
        
        #print("+output_dir="+output_dir)

        cmdline = " ".join(arguments) #join all arguments passed that aren't number of tasks or task name
        # installation.yaml's own prosculpt_python_path, not sys.executable
        # and not a bare "python": sys.executable is whatever interpreter
        # happens to be running *this* script, which depends on how the
        # caller invoked slurm_runner.py - fine when the dashboard does it
        # (it explicitly launches slurm_runner.py with its configured
        # python_path), but silently wrong if a human runs `python
        # slurm_runner.py ...` without having activated the prosculpt conda
        # env first (or any other caller whose active interpreter isn't
        # prosculpt's). A bare "python" has the matching problem one step
        # later: it gets resolved fresh on the compute node when the job
        # actually runs, which isn't guaranteed to have the right env either.
        # installation.yaml's prosculpt_python_path is the one absolute path
        # meant to always be correct regardless of who/what invoked this
        # script - it's the same value prosculpt_run.py itself uses to call
        # scoring_script.py, so this keeps both call sites in sync.
        line = f"""{installation_yaml_data['prosculpt_python_path']} {slurm_runner_path}/prosculpt_run.py {cmdline}"""
        print(line, file=f)

print(f"Slurm command can be found in {out_command_file}")

options_string=""
job_specific_keys=[]
if slurm_data is not None:
    for key, value in slurm_data.items(): 
        job_specific_keys.append(key)
        # Recorded in job_specific_keys either way (so the installation.yaml
        # loop below doesn't also try to add its own output/error), but not
        # actually emitted here when forced_log_dir is set - that's the
        # whole point of --allow-custom-log-path defaulting to off: the
        # job yaml's own output/error choice is exactly what gets
        # overridden by default.
        if key in ("output", "error") and forced_log_dir is not None:
            continue
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
        if key in ("output", "error") and forced_log_dir is not None:
            continue
        if len(key)==1:
            options_string+= " -"
        else:
            options_string+= " --"
        options_string+= f"{key} {value}"
        

# The fallback here (relative "logs/slurm-%A_%a.err"/".out", same as
# always) only actually applies with --allow-custom-log-path - otherwise
# it's immediately overridden by whichever of the job yaml's or
# installation.yaml's own output/error survived the filtering above (both
# come later in the command line than these two flags, and the last -o/-e
# sbatch sees wins), same as before this option existed at all. With
# forced_log_dir set, that filtering already removed both of those, so
# these are the only -o/-e sbatch actually sees.
log_err = f"{forced_log_dir}/slurm-%A_%a_%x.err" if forced_log_dir is not None else "logs/slurm-%A_%a.err"
log_out = f"{forced_log_dir}/slurm-%A_%a_%x.out" if forced_log_dir is not None else "logs/slurm-%A_%a.out"

if not args[0].dry_run:
    # --parsable makes sbatch print just the numeric job id (plus ";cluster"
    # on a federated setup) to stdout instead of its normal "Submitted batch
    # job N" banner - captured here (rather than the os.system() this used
    # to be, which discards output entirely) so callers that need the real
    # job id (e.g. the dashboard's job-submission API, to know which log
    # file to watch for) can read it reliably instead of scraping free-form
    # text.
    full_command= f"export GROUP_SIZE=1; sbatch --parsable -J {task_name} -a 1-{n} -e {log_err} -o {log_out}  {options_string} {slurm_runner_path}/wrapper_slurm_array_job_group.sh {out_command_file}"
    print(f"Full command is: {full_command}")
    result = subprocess.run(full_command, shell=True, capture_output=True, text=True)
    exit_code = result.returncode
    if result.stdout:
        print(result.stdout, end="" if result.stdout.endswith("\n") else "\n")
    if result.stderr:
        print(result.stderr, end="" if result.stderr.endswith("\n") else "\n", file=sys.stderr)
    if exit_code == 0:
        job_id = result.stdout.strip().split(";")[0]
        print(f"Job {task_name} has been submitted to slurm with id {job_id} and code {exit_code}")
    else:
        print(f"Job submission failed with code {exit_code}")
    # Propagate sbatch's own exit code as this script's exit code - without
    # this, the block above already prints "Job submission failed..." on a
    # real failure (e.g. sbatch couldn't reach the controller) but the
    # script itself still finishes and exits 0 regardless, which makes a
    # genuine submission failure look identical to success to anything
    # checking $? - a plain shell caller, and in particular the dashboard's
    # own submit API, which reports "Submitted" (ok: true) purely from the
    # remote command's exit code.
    sys.exit(exit_code)
else:
    print("Command wasn't run because --dry-run was active.")
