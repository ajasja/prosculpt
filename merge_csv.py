
import os
import pandas as pd
import argparse
import help_functions

parser = argparse.ArgumentParser(description="Run protein design pipeline")
parser.add_argument("--output_dir", help="Path to output directory")

args = parser.parse_args()



def merge_csv_files(output_dir):
    
    merged_data_frames = [] #cannot merge df directly because df.append is depreciated

    for sub_directory in os.listdir(output_dir): #os.listdir() returns a list of all the files and directories in a specified directory.
        subdir_path = os.path.join(output_dir, sub_directory)
        if os.path.isdir(subdir_path): #if an element of the list is a directory:
            for filename in os.listdir(subdir_path):
                if filename.endswith(".csv"):
                    file_path = os.path.join(subdir_path, filename)
                    data = pd.read_csv(file_path)
                    merged_data_frames.append(data)

    merged_data = pd.concat(merged_data_frames, ignore_index=True)
    merged_file_path = os.path.join(output_dir, "merged.csv")
    merged_data.to_csv(merged_file_path, index=False)
    print(f"Merged CSV file saved at: {merged_file_path}")


merge_csv_files(args.output_dir)

