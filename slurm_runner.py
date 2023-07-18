import sys
import shlex
import time
import os

WIDTH = 80
SYMBOL = "⁝"
print("".center(WIDTH, SYMBOL))
print(" prosculpt2slurm ".center(WIDTH, SYMBOL))
print("".center(WIDTH, SYMBOL))
print("Please make sure your first argument was the number of concurrent jobs you want to run and the second argument was the desired name of the slurm task. After that, add the arguments that will be passed as-is to the _merged script.")
print("".center(WIDTH, SYMBOL))


n = int(sys.argv[1])
task_name = sys.argv[2]

if task_name[:11] == "output_dir=":
    print(f"!!! Are you sure you want to name your slurm job {task_name}? Anyway, continuing ...")

out_command_file = f"ps2slurm_{int(time.time())}.txt"

with open(out_command_file, 'w') as f:
    for i in range(n):
        let_argsv = sys.argv.copy()
        let_argsv[3] += f"/{i:02d}"

        cmdline = " ".join(map(shlex.quote, let_argsv[3:]))
        line = f"""python /home/zznidar/dn/prs/NOVO/rfdiff_mpnn_af2_merged.py {cmdline}"""
        print(line, file=f)

print(f"Slurm command can be found in {out_command_file}")
os.system(f"export GROUP_SIZE=1; sbatch --partition=gpu --gres=gpu:A40:1 --ntasks=1 --cpus-per-task=2 -J {task_name}  -a 1-{n} wrapper_slurm_array_job_group.sh {out_command_file}")
print(f"Job {task_name} has been submitted to slurm".center(WIDTH, SYMBOL))