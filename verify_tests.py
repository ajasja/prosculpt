import argparse
import os
import pandas as pd 
import datetime

parser=argparse.ArgumentParser()
parser.add_argument("dir", help='Directory of test outputs.')
args=parser.parse_args()

directory=args.dir

with open(os.path.join(directory,"test_verification_output.txt"), "a+") as output_file:
    for test_folder in os.listdir(directory):
        if os.path.isdir(os.path.join(directory,test_folder)):
            final_output_path=os.path.join(directory,test_folder,"00","final_output.csv")
            final_pdbs_path=os.path.join(directory,test_folder,"final_pdbs")

            if os.path.exists(final_output_path):
                final_output=pd.read_csv(final_output_path)
                if len(final_output)==len(os.listdir(final_pdbs_path)):
                    output_file.write("Test " + test_folder + " passed.\n")
                else:
                    output_file.write("Test " + test_folder + "failed. final_output.csv has "+str(len(final_output))+ " lines but only "+
                        len(os.listdir(final_pdbs_path)) + " files in final_pdbs folder\n")
            else:
                output_file.write("Test "+ test_folder + " failed. Missing final_output.csv\n")

    output_file.write("Finished running and verifying tests: "+datetime.datetime.now().strftime('%d/%m/%y %H:%M:%S.%f')+"\n")
