

n = 2

f"""/home/nbizjak/projects/11_04_2023_rigid_connections/rfdiff_mpnn_af2_disontiuous.py --pdb_path /home/nbizjak/projects/11_04_2023_rigid_connections/inputs_for_rfdiff/2442_ver1_chains-trimm.pdb --output_dir /home/nbizjak/projects/29_05_2023_zvezda_test/outputs/ --contig [E10-68/5-15/C4-72/5-15/D78-145/5-15/A1-60/10/B63-119/5-15/F86-144] --num_designs_rfdiff 2 --num_seq_per_target_mpnn 5 --chains_to_design A"""

output_file = 'for_slurm_prosculpt_dev.txt'


with open(output_file, 'w') as f:
    for i in range(n):
        
        line = f"""/home/nbizjak/projects/11_04_2023_rigid_connections/.venv/bin/python /home/nbizjak/prosculpt_dev/prosculpt/rfdiff_mpnn_af2_disontiuous.py --original_pdb_path "/home/nbizjak/projects/11_04_2023_rigid_connections/inputs_for_rfdiff/mALb8-antiA.pdb" --output_dir  /home/nbizjak/prosculpt_dev/test_ouptuts_dir/outputs/{i:02d} --contig '[C33-60/5-7/A1-30/0 B61-120]' --num_designs_rfdiff 2  --num_seq_per_target_mpnn 2 --chains_to_design_mpnn 'A' --af2_mpnn_cycles 2"""
        print(line, file=f)

"""
         line = params["path_to_script"]
        for key, value in params.items():
            if not key=="path_to_script":
                line += f" --{key} {value}"
"""
# os system
##export GROUP_SIZE=1; sbatch --partition=gpu --gres=gpu:A40:1 --ntasks=1 --cpus-per-task=2 -J zvezda  -a 1-20 /home/aljubetic/scripts/wrapper_slurm_array_job_group.sh 35_06_2023_brez_simetrije_run1/for_slurm_brez_sim_run1.txt


#export GROUP_SIZE=1; sbatch --partition=gpu --exclude=compute-6-0,compute-6-2,compute-6-4,compute-3-19,compute-3-20  --gres=gpu:A40:1 --ntasks=1 --cpus-per-task=2 -J zvezda  -a 1-60%20 /home/aljubetic/scripts/wrapper_slurm_array_job_group.sh for_slurm_zvezda.txt

#export GROUP_SIZE=1; sbatch --partition=gpu --exclude=compute-6-0,compute-6-2,compute-6-4,compute-3-19,compute-3-20  --gres=gpu:A40:1 --ntasks=1 --cpus-per-task=2 -J prosculpt_dev  -a 1-20 /home/aljubetic/scripts/wrapper_slurm_array_job_group.sh for_slurm_prosculpt_dev.txt