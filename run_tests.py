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
os.makedirs("Examples/Examples_out/", exist_ok=True)

if args.list:
    print("Available tests:")
    for sh_file in glob.glob("Examples/*.sh"):
        print(sh_file.split("/")[1].split(".")[0])   
    sys.exit()

#Assign either one test or all to the list
test_file_list=[]
if args.test:
    test_filename="Examples/"+args.test+".sh"
    assert os.path.exists(test_filename), "Test not found. Use the name without the extension."
    test_file_list.append(test_filename)
else:
    for sh_file in glob.glob("Examples/*.sh"):
        test_file_list.append(sh_file)


#then run the list of tests
slurm_job_list=[]
for test_file in test_file_list:
    print("Running "+ test_file)
    if not args.dry_run:
        process_output=subprocess.run(["sh",test_file],capture_output=True,text=True)
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
    os.system(f"sbatch -d afterany{jobs_string} -J test_check slurm_verify_tests.sh")
    

with open("Examples/Examples_out/test_verification_output.txt","w+") as output_file:
    output_file.write("Started running tests: "+datetime.datetime.now().strftime('%d/%m/%y %H:%M:%S.%f')+"\n")

            


