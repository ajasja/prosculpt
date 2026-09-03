import argparse
import os
import re
import pandas as pd
import yaml
import datetime

print("Running test verification script: "+datetime.datetime.now().strftime('%d/%m/%y %H:%M:%S.%f'))

parser=argparse.ArgumentParser()
parser.add_argument("dir", help='Directory of test outputs.')
args=parser.parse_args()

directory=args.dir

print(f"Verifying tests in {directory}")

# Subfolder names prosculpt creates inside a task's output_dir, not a test folder.
RESERVED_DIR_NAMES = {"output_pdbs", "filtered_pdbs", "logs", "scoring_logs"}

# A SLURM array task's own numbered subdirectory (not fixed-width - task 100 is "100").
_TASK_DIR_RE = re.compile(r"^\d+$")


def find_task_dirs(test_dir):
    """Every numbered task subdirectory of a test's output directory, sorted
    numerically. Falls back to the test directory itself if there are none."""
    task_dirs = []
    for name in os.listdir(test_dir):
        if _TASK_DIR_RE.match(name) and os.path.isdir(os.path.join(test_dir, name)):
            task_dirs.append((int(name), os.path.join(test_dir, name)))
    task_dirs.sort(key=lambda t: t[0])
    if task_dirs:
        return [d for _, d in task_dirs]
    return [test_dir]


def post_filtering_scoring_configured(task_dir):
    """Whether this task's input.yaml requested a post_filtering_scoring_script."""
    input_yaml_path = os.path.join(task_dir, "input.yaml")
    if not os.path.isfile(input_yaml_path):
        return False
    try:
        with open(input_yaml_path) as f:
            data = yaml.safe_load(f) or {}
    except yaml.YAMLError:
        return False
    filtering_cfg = data.get("filtering") or {}
    return bool(filtering_cfg.get("post_filtering_scoring_script"))


def verify_task(task_dir):
    """Returns (passed, message) for one task directory: checks output.csv
    exists, filtered_output.csv exists, and (if post-filtering scoring was
    configured) that filtered_output.csv gained columns beyond output.csv."""
    output_csv_path = os.path.join(task_dir, "output.csv")
    if not os.path.isfile(output_csv_path):
        return False, "Missing output.csv"

    filtered_csv_path = os.path.join(task_dir, "filtered_output.csv")
    if not os.path.isfile(filtered_csv_path):
        return False, "Missing filtered_output.csv (filtering stage did not run/complete)"

    if post_filtering_scoring_configured(task_dir):
        try:
            # nrows=0: just the header.
            output_columns = set(pd.read_csv(output_csv_path, nrows=0).columns)
            filtered_columns = set(pd.read_csv(filtered_csv_path, nrows=0).columns)
        except Exception as e:
            return False, f"Could not read output.csv/filtered_output.csv to check post-filtering scoring: {e}"
        if not (filtered_columns - output_columns):
            return False, (
                "post_filtering_scoring_script was configured, but filtered_output.csv has no "
                "columns beyond output.csv - the post-filtering scoring job may not have completed"
            )

    return True, "passed"


with open(os.path.join(directory,"test_verification_output.txt"), "a+") as output_file:
    print("Opened test_verification_output.txt")
    for test_folder in sorted(os.listdir(directory)):
        if test_folder in RESERVED_DIR_NAMES:
            continue
        test_dir = os.path.join(directory, test_folder)
        if not os.path.isdir(test_dir):
            continue
        print(f"test_folder: {test_folder}")

        for task_dir in find_task_dirs(test_dir):
            # Label each replica distinctly, e.g. "binder_AF3/01".
            task_label = test_folder if task_dir == test_dir else f"{test_folder}/{os.path.basename(task_dir)}"
            passed, message = verify_task(task_dir)
            if passed:
                print(f"{task_label}: passed")
                output_file.write("Test " + task_label + " passed.\n")
            else:
                print(f"{task_label}: failed - {message}")
                output_file.write("Test " + task_label + " failed. " + message + "\n")

    output_file.write("Finished running and verifying tests: "+datetime.datetime.now().strftime('%d/%m/%y %H:%M:%S.%f')+"\n")
print("Finished running test verification script: "+datetime.datetime.now().strftime('%d/%m/%y %H:%M:%S.%f'))