

n = 20

f"""/home/nbizjak/projects/11_04_2023_rigid_connections/rfdiff_mpnn_af2_disontiuous.py --pdb_path /home/nbizjak/projects/11_04_2023_rigid_connections/inputs_for_rfdiff/2442_ver1_chains-trimm.pdb --output_dir /home/nbizjak/projects/29_05_2023_zvezda_test/outputs/ --contig [E10-68/5-15/C4-72/5-15/D78-145/5-15/A1-60/10/B63-119/5-15/F86-144] --num_designs_rfdiff 2 --num_seq_per_target_mpnn 5 --chains_to_design A"""

output_file = 'for_slurm_C11_sim_.txt'


with open(output_file, 'w') as f:
    for i in range(n):
        
        line = f"""/home/nbizjak/projects/11_04_2023_rigid_connections/.venv/bin/python /home/nbizjak/projects/30_05_2023_simetricni_motv_test/rfdiff_mpnn_af2_simetrija.py --original_pdb_path /home/nbizjak/projects/30_05_2023_simetricni_motv_test/1qaw-min-aligned.pdb --output_dir  /home/nbizjak/projects/40_05_2023_C11_simetrija_run_1/outputs/{i:02d} --contig '[15/A29-41/15/0 15/B29-41/15/0 15/C29-41/15/0 15/D29-41/15/0 15/E29-41/15/0 15/F29-41/15/0 15/G29-41/15/0 15/H29-41/15/0 15/I29-41/15/0 15/J29-41/15/0 15/K29-41/15]' --num_designs_rfdiff 15  --num_seq_per_target_mpnn 10 --chains_to_design_mpnn 'A B C D E F G H I J K' """
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

#export GROUP_SIZE=1; sbatch --partition=gpu --exclude=compute-6-0,compute-6-2,compute-6-4,compute-3-19,compute-3-20  --gres=gpu:A40:1 --ntasks=1 --cpus-per-task=2 -J C11_sim_run1  -a 1-20 /home/aljubetic/scripts/wrapper_slurm_array_job_group.sh for_slurm_C11_sim_.txt