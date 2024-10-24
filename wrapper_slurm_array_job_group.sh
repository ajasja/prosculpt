#!/bin/bash
#SBATCH -p intel 
#SBATCH --mem=4g 

GROUP_SIZE=${GROUP_SIZE:-1}
echo "GROUP SIZE IS: $GROUP_SIZE"

export NVIDIA_VISIBLE_DEVICES=$GPU_DEVICE_ORDINAL
export CUDA_VISIBLE_DEVICES=$GPU_DEVICE_ORDINAL

for I in $(seq 1 $GROUP_SIZE)
do
    echo "Hello from job $SLURM_JOB_ID on $(hostname) at $(date)"
    echo "PWD:"
    pwd
    echo "which python? $(which python)"
    echo ""
    J=$(($SLURM_ARRAY_TASK_ID * $GROUP_SIZE + $I - $GROUP_SIZE))
    CMD=$(sed -n "${J}p" $1)
    echo "COMMAND: ${CMD}"
    echo "${CMD}" | bash
done


#####command:
#####export GROUP_SIZE={GROUP_SIZE}; sbatch --mem=4G -p short -J {task_list} -o {task_list}.out -e {task_list}.err -a 1-{num_tasks} /home/aljubetic/scripts/wrapper_slurm_array_job_group.sh {task_list}




