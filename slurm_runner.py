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
        "By default, sbatch's own -o/-e are forced to <output_dir>/logs/ (and "
        "<output_dir>/scoring_logs/ for the post-filtering scoring job), overriding any "
        "slurm.output/error or filtering_slurm.output/error set in the job yaml or "
        "installation.yaml. Pass this to use those settings instead."
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

# filtering: is optional. post_filtering_scoring_script enables a second,
# dependent SLURM job below; filtering_slurm holds that job's slurm options,
# merged against installation.yaml's filtering_slurm defaults.
filtering_data = yaml_data.get("filtering") or {}
post_filtering_scoring_script = filtering_data.get("post_filtering_scoring_script")
filtering_slurm_data = filtering_data.get("filtering_slurm") or {}
installation_filtering_slurm_yaml_data = installation_yaml_data.get("filtering_slurm") or {}

if "task_name" not in yaml_data:
    task_name="prosculpt_task"
else:
    task_name=yaml_data["task_name"]

# Base output_dir - what the per-task loop below builds each task's numbered
# dir from. Also used as where sbatch's forced -o/-e and the generated
# ps2slurm_*.txt command files go. Must stay the *base* dir, not a per-task
# one: all tasks share a single sbatch array submission, and SLURM's %a
# substitution isn't zero-padded to match "01"/"02"/..., so a shared
# template can't target each task's own numbered directory.
base_output_dir = None
for arg in extra_args:
    if "output_dir" in arg:
        base_output_dir = arg.split("=", 1)[-1]
if base_output_dir is None:
    base_output_dir = yaml_data["output_dir"]

# Needed even with --allow-custom-log-path, since the command file below is
# written directly into this directory.
os.makedirs(base_output_dir, exist_ok=True)

out_command_file = os.path.join(base_output_dir, f"ps2slurm_{task_name}_{int(time.time_ns())}.txt")


def resolve_forced_log_dir(subdir_name):
    """
    Returns base_output_dir/subdir_name (creating it if needed), or None if
    --allow-custom-log-path was passed.
    """
    if args[0].allow_custom_log_path:
        return None
    forced_dir = os.path.join(base_output_dir, subdir_name)
    # sbatch does not create the directory portion of -o/-e itself - if it
    # doesn't already exist when the job is submitted, the job fails
    # immediately with no log written anywhere at all (see
    # job_staging.ensure_logs_dir() in the dashboard for the same problem
    # solved the same way, for the cwd-relative logs/ this replaces).
    os.makedirs(forced_dir, exist_ok=True)
    return forced_dir


forced_log_dir = resolve_forced_log_dir("logs")


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


def build_options_string(job_specific_data, installation_default_data, forced_log_dir):
    """
    Builds sbatch's extra options string from a job-specific slurm-like dict
    merged against installation.yaml defaults of the same shape; job-specific
    keys win. If forced_log_dir is set, "output"/"error" keys are dropped
    from both (the caller supplies its own -o/-e instead).
    """
    options_string = ""
    job_specific_keys = []
    for key, value in (job_specific_data or {}).items():
        job_specific_keys.append(key)
        if key in ("output", "error") and forced_log_dir is not None:
            continue
        if key == "slurm_options_string":
            options_string += f" {value}"
        else:
            options_string += " -" if len(key) == 1 else " --"
            options_string += f"{key} {value}"

    for key, value in (installation_default_data or {}).items():
        if key in job_specific_keys:
            continue
        if key in ("output", "error") and forced_log_dir is not None:
            continue
        options_string += " -" if len(key) == 1 else " --"
        options_string += f"{key} {value}"

    return options_string


options_string = build_options_string(slurm_data, installation_slurm_yaml_data, forced_log_dir)

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

exit_code = 0
main_job_id = None
if not args[0].dry_run:
    # --parsable makes sbatch print just the job id instead of its normal
    # banner, so callers (the dashboard, and the post-filtering scoring
    # job's --dependency=aftercorr:<main_job_id> below) can read it reliably.
    full_command= f"export GROUP_SIZE=1; sbatch --parsable -J {task_name} -a 1-{n} -e {log_err} -o {log_out}  {options_string} {slurm_runner_path}/wrapper_slurm_array_job_group.sh {out_command_file}"
    print(f"Full command is: {full_command}")
    result = subprocess.run(full_command, shell=True, capture_output=True, text=True)
    exit_code = result.returncode
    if result.stdout:
        print(result.stdout, end="" if result.stdout.endswith("\n") else "\n")
    if result.stderr:
        print(result.stderr, end="" if result.stderr.endswith("\n") else "\n", file=sys.stderr)
    if exit_code == 0:
        main_job_id = result.stdout.strip().split(";")[0]
        print(f"Job {task_name} has been submitted to slurm with id {main_job_id} and code {exit_code}")
    else:
        print(f"Job submission failed with code {exit_code}")
else:
    print("Command wasn't run because --dry-run was active.")

# Post-filtering scoring: a second array job, submitted only if
# post_filtering_scoring_script is set. Depends on the main job via
# --dependency=aftercorr:<main_job_id> - each array index starts once that
# same index in the main job completes, independently of the rest.
if post_filtering_scoring_script:
    scoring_command_file = os.path.join(base_output_dir, f"ps2slurm_{task_name}_scoring_{int(time.time_ns())}.txt")
    with open(scoring_command_file, 'w') as f:
        for i in range(1, n + 1):
            arguments = []

            output_dir = ""
            output_dir_in_args = False
            for arg in extra_args:
                if "output_dir" in arg:
                    output_dir_in_args = True
                    output_dir = arg
                    if output_dir[-1:] != "/":
                        output_dir += "/"
                    output_dir += f"{i:02d}"
                    arguments.append(output_dir)

            if not output_dir_in_args:
                output_dir = yaml_data["output_dir"]
                if output_dir[-1:] != "/":
                    output_dir += "/"
                output_dir += f"{i:02d}"
                arguments.append("++output_dir=" + output_dir)

            # Must come before -cd/-cn, like any other hydra override.
            arguments.append("+only_run_post_filtering_scoring=true")
            arguments.append(f"-cd '{yaml_file_dir}'")
            arguments.append(f"-cn '{yaml_file_name_no_extension}'")

            cmdline = " ".join(arguments)
            line = f"""{installation_yaml_data['prosculpt_python_path']} {slurm_runner_path}/prosculpt_run.py {cmdline}"""
            print(line, file=f)

    print(f"Post-filtering scoring slurm command can be found in {scoring_command_file}")

    forced_scoring_log_dir = resolve_forced_log_dir("scoring_logs")
    scoring_options_string = build_options_string(
        filtering_slurm_data, installation_filtering_slurm_yaml_data, forced_scoring_log_dir
    )
    scoring_log_err = (
        f"{forced_scoring_log_dir}/slurm-%A_%a_%x.err"
        if forced_scoring_log_dir is not None
        else "scoring_logs/slurm-%A_%a.err"
    )
    scoring_log_out = (
        f"{forced_scoring_log_dir}/slurm-%A_%a_%x.out"
        if forced_scoring_log_dir is not None
        else "scoring_logs/slurm-%A_%a.out"
    )

    if args[0].dry_run:
        placeholder_dependency = "--dependency=aftercorr:<main_job_id>  # not known in a dry run - only real once the main job is actually submitted"
        scoring_full_command = f"sbatch --parsable -J {task_name}_scoring -a 1-{n} {placeholder_dependency} -e {scoring_log_err} -o {scoring_log_out} {scoring_options_string} {slurm_runner_path}/wrapper_slurm_array_job_group.sh {scoring_command_file}"
        print(f"Post-filtering scoring command (dry run): {scoring_full_command}")
        print("Post-filtering scoring command wasn't run because --dry-run was active.")
    elif main_job_id is None:
        print(
            "Skipping post-filtering scoring submission: the main job was not submitted successfully."
        )
    else:
        scoring_full_command = f"sbatch --parsable -J {task_name}_scoring -a 1-{n} --dependency=aftercorr:{main_job_id} -e {scoring_log_err} -o {scoring_log_out} {scoring_options_string} {slurm_runner_path}/wrapper_slurm_array_job_group.sh {scoring_command_file}"
        print(f"Full post-filtering scoring command is: {scoring_full_command}")
        scoring_result = subprocess.run(scoring_full_command, shell=True, capture_output=True, text=True)
        scoring_exit_code = scoring_result.returncode
        if scoring_result.stdout:
            print(scoring_result.stdout, end="" if scoring_result.stdout.endswith("\n") else "\n")
        if scoring_result.stderr:
            print(scoring_result.stderr, end="" if scoring_result.stderr.endswith("\n") else "\n", file=sys.stderr)
        if scoring_exit_code == 0:
            scoring_job_id = scoring_result.stdout.strip().split(";")[0]
            print(
                f"Post-filtering scoring job {task_name}_scoring has been submitted to slurm with id "
                f"{scoring_job_id} and code {scoring_exit_code}, depending on job {main_job_id}"
            )
        else:
            print(f"Post-filtering scoring job submission failed with code {scoring_exit_code}")
        if scoring_exit_code != 0:
            exit_code = scoring_exit_code

# Exit with the first non-zero submission's exit code (main or scoring), so
# a real submission failure isn't reported as success to callers checking $?.
if not args[0].dry_run:
    sys.exit(exit_code)
