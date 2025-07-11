import os
import argparse
import glob
import sys
import subprocess
import datetime

parser=argparse.ArgumentParser()

parser.add_argument('-t', '--test', help='Name of test to run (without extension). If not specified, all tests will be ran.')
parser.add_argument('-ls', '--list', action="store_true", help="See list of available tests.")
parser.add_argument('-d', '--dry-run', action="store_true", help="Do not run anything. Mostly for debugging")


args=parser.parse_args()
out_folder=f"Examples_out_{datetime.datetime.now().strftime('%y-%m-%d-%H-%M-%S')}"
os.makedirs(f"Examples/{out_folder}", exist_ok=True)

if args.list:
    print("Available tests:")
    for sh_file in glob.glob("Examples/*.yaml"):
        print(sh_file.split("/")[1].split(".")[0])   
    sys.exit()

#Assign either one test or all to the list
test_file_list=[]
if args.test:
    test_filename="Examples/"+args.test+".yaml"
    assert os.path.exists(test_filename), "Test not found. Use the name without the extension."
    test_file_list.append(test_filename)
else:
    for sh_file in glob.glob("Examples/*.yaml"):
        test_file_list.append(sh_file)


#then run the list of tests
slurm_job_list=[]
for test_file in test_file_list:
    print("Running "+ test_file)
    if not args.dry_run:
        
        command=f"python slurm_runner.py {test_file} ++output_dir='Examples/{out_folder}/{test_file.split("/")[1].split(".")[0]}'"
        print(command)
        process_output=subprocess.run(command, shell=True,capture_output=True,text=True)
        for line in process_output.stdout.split("\n"):
            if "Submitted batch job" in line:
                print(line)
                slurm_job_list.append(line.split("job ")[1])
        if process_output.returncode!=0:
            print(f"ERROR {test_file}:")
            for line in process_output.stderr.split("\n"):
                print(line)

print(f"Slurm job list: {slurm_job_list}")
#Queue up the test verifier after all jobs have ended
jobs_string=""
for id, job in enumerate(slurm_job_list):
    jobs_string+=":"
    jobs_string+=job
        
#subprocess.run(command)
if not args.dry_run:
    os.system(f"sbatch -d afterany{jobs_string} -J test_check slurm_verify_tests.sh Examples/{out_folder}") # compute-0-5 fails to run python (command not found, although PATH is correct). compute-0-2 and compute-0-10 also fail (python not found; if I provide the full path to my .venv python, no such file or directory). Thus, we exclude all compute-0- nodes from 1 to 10.
    with open(f"Examples/{out_folder}/test_verification_output.txt","w+") as output_file:
        output_file.write("Started running tests: "+datetime.datetime.now().strftime('%d/%m/%y %H:%M:%S.%f')+"\n")

                


