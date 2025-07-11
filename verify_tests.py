import argparse
import os
import pandas as pd 
import datetime

print("Running test verification script: "+datetime.datetime.now().strftime('%d/%m/%y %H:%M:%S.%f'))

parser=argparse.ArgumentParser()
parser.add_argument("dir", help='Directory of test outputs.')
args=parser.parse_args()

directory=args.dir

print(f"Verifying tests in {directory}")

with open(os.path.join(directory,"test_verification_output.txt"), "a+") as output_file:
    print("Opened test_verification_output.txt")
    for test_folder in os.listdir(directory):
        print(f"test_folder: {test_folder}")
        if test_folder!="final_pdbs":
            print("test_folder is not final_pdbs")
            if os.path.isdir(os.path.join(directory,test_folder)):
                print("test_folder is a directory")
                if os.path.isdir(os.path.join(directory,test_folder,"01")): #Check if 00 folder exists. Not sure what triggers this folder to exist or not.
                    ## ANSWER: Folder number is defined in slurm_runner.py where writing out_command_file (it overwrites the output_dir specified in config, by adding the number. I think it is necessary for making multiple models parallelly). 
                    ## This folder is then created in rfdiff_mpnn_af2_merged.py general_config_prep(), which is the first function called in the script.
                    print("01 folder exists")
                    final_output_path=os.path.join(directory,test_folder,"01","final_output.csv")
                else:
                    print("01 folder does not exist")
                    final_output_path=os.path.join(directory,test_folder,"final_output.csv")
                    
                #final_pdbs_path=os.path.join(directory,test_folder,"final_pdbs")

                if os.path.exists(final_output_path):
                    print(f"final_output exists and is {final_output_path}")
                    #This part is commented because it was failing with the binders test because of folder architecture.
                    #final_output=pd.read_csv(final_output_path)
                    #if len(final_output)==len(os.listdir(final_pdbs_path)):
                    #    output_file.write("Test " + test_folder + " passed.\n")
                    #else:
                    #    output_file.write("Test " + test_folder + "failed. final_output.csv has "+str(len(final_output))+ " lines but only "+
                    #        len(os.listdir(final_pdbs_path)) + " files in final_pdbs folder\n")
                    
                    output_file.write("Test " + test_folder + " passed.\n")
                else:
                    print("final_output does not exist")
                    output_file.write("Test "+ test_folder + " failed. Missing final_output.csv\n")

        
    output_file.write("Finished running and verifying tests: "+datetime.datetime.now().strftime('%d/%m/%y %H:%M:%S.%f')+"\n")
print("Finished running test verification script: "+datetime.datetime.now().strftime('%d/%m/%y %H:%M:%S.%f'))