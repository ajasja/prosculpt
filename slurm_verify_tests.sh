#!/bin/sh
echo "Hello from job $SLURM_JOB_ID on $(hostname) at $(date). I am verifying the tests."
python verify_tests.py Examples/Examples_out/