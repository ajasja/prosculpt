#!/bin/sh
echo "Hello from job $SLURM_JOB_ID on $(hostname) at $(date). I am verifying the tests."
echo $PATH 
python verify_tests.py $1