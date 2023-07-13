print("Please note that Hydra, OmegaConf etc. are now required.")

n = 2


output_file = 'ab_cdr_A.txt'


with open(output_file, 'w') as f:
    for i in range(n):
        
        #line = f"""/home/nbizjak/projects/11_04_2023_rigid_connections/.venv/bin/python /home/zznidar/dn/prs/rfdiff_mpnn_af2_disontiuous.py --output_dir  ZZtest/ab_cdr/STAR_out_TEST/{i:02d} --contig \[H1-25/8-8/H34-51/6-6/H58-98/22-22/H121-127\] --original_pdb_path "/home/zznidar/dn/prs/tetetetettest.pdb" --num_designs_rfdiff 1  --num_seq_per_target_mpnn 1 --chains_to_design_mpnn H --af2_mpnn_cycles 1 """

        #NOVO: ### NOTE: It requires you to have python installed and be in .venv activated with all modules installed.
        line = f"""python /home/zznidar/dn/prs/NOVO/rfdiff_mpnn_af2_merged.py output_dir=ZZtest/ab_cdr/NOVO_out_9/{i:02d} "contig=\[H1-25/8-8/H34-51/6-6/H58-98/22-22/H121-127\]" pdb_path="/home/zznidar/dn/prs/tetetetettest.pdb" num_designs_rfdiff=2 num_seq_per_target_mpnn=2 chains_to_design=A af2_mpnn_cycles=2"""
        
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

#export GROUP_SIZE=1; sbatch --partition=gpu --exclude=compute-6-0,compute-6-2,compute-6-4,compute-3-19,compute-3-20  --gres=gpu:A40:1 --ntasks=1 --cpus-per-task=2 -J prosculpt_dev  -a 1-2 /home/aljubetic/scripts/wrapper_slurm_array_job_group.sh for_slurm_two_helices_3cycles_model4.txt






# python /home/zznidar/dn/prs/NOVO/rfdiff_mpnn_af2_merged.py output_dir=ZZtest/ab_cdr/NOVO_out_1/i contig="\[H1-25/8-8/H34-51/6-6/H58-98/22-22/H121-127\]" pdb_path="/home/zznidar/dn/prs/tetetetettest.pdb" num_designs_rfdiff=1 num_seq_per_target_mpnn=1 chains_to_design='H' af2_mpnn_cycles=1


"""
export GROUP_SIZE=1; sbatch --partition=gpu --gres=gpu:A40:1 --ntasks=1 --cpus-per-task=2 -J name_of_task  -a 1-2 wrapper_slurm_array_job_group.sh ab_cdr_deb.txt
"""