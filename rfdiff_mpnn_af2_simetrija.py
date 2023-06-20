import os
import help_functions
import glob
import argparse
import re
import shutil

"""
READ ME
Guidelines:
    - Script written to design MULTIPLE linkers INSIDE ONE CHAIN or ONE linker to connect TWO CHAINS
    - RFdiff designs protein so the designed chain/chains are now chain A, regardless of how they were identified before
        - Other undesigned chains are B, C,... regardles if before they were A

"""

#/home/nbizjak/projects/11_04_2023_rigid_connections/inputs_for_rfdiff/phosphoCC_bundle.pdb
#/home/nbizjak/projects/11_04_2023_rigid_connections/inputs_for_rfdiff/dva_helixa_delujoca.pdb
pdb_path ="/home/nbizjak/projects/30_05_2023_simetricni_motv_test/1qaw-min-aligned.pdb"
output_dir = "/home/nbizjak/projects/30_05_2023_simetricni_motv_test/outputs/"
contig = "[10/A29-41/10/0 10/B29-41/10/0 10/C29-41/10/0 10/D29-41/10/0 10/E29-41/10/0 10/F29-41/10/0 10/G29-41/10/0 10/H29-41/10/0 10/I29-41/10/0 10/J29-41/10/0 10/K29-41/10]" #[A1-30/5/C1-30/0 D1-30/0  B1-30] #[A1-37/30-60/A42-70] dva helixa  #[C33-60/3/A1-30/0 B61-120] mlba... [C33-60/5/A1-30/0 B61-120] ...[10/A29-41/10/0 10/B29-41/10/0 10/C29-41/10/0 10/D29-41/10/0 10/E29-41/10/0 10/F29-41/10/0 10/G29-41/10/0 10/H29-41/10/0 10/I29-41/10/0 10/J29-41/10/0 10/K29-41/10]
num_designs_rfdiff = 1
num_seq_per_target_mpnn = 1
chains_to_design = "A B C D E F G H I J K" #might be very important for complexes


#if designing complexes change in helper funcitions in the pdb proceisng function

parser = argparse.ArgumentParser(description="Run protein design pipeline")
parser.add_argument("--contig", default=contig, help="Contig to design")
parser.add_argument("--original_pdb_path", default=pdb_path, help="Path to input PDB file")
parser.add_argument("--output_dir", default = output_dir, help="Path to output directory")
parser.add_argument("--num_designs_rfdiff", type=int, default=num_designs_rfdiff, help="Number of designs to generate with RFdiffusion")
parser.add_argument("--num_seq_per_target_mpnn", type=int, default=num_seq_per_target_mpnn, help="Number of sequences to generate per target with MPNN")
parser.add_argument("--chains_to_design_mpnn", default= chains_to_design, help="All chains of protein even if not redesigned (except if doing binder design)") #could be infered from pdb

args = parser.parse_args()

#_____________ GENERAL PREP ______________

os.makedirs(args.output_dir, exist_ok=True) #exist_ok=True ensures the function does not raise an error if the directory already exists
rfdiff_out_dir = os.path.join(args.output_dir, "1_rfdiff")
mpnn_out_dir = os.path.join(args.output_dir, "2_mpnn")
af2_out_dir = os.path.join(args.output_dir, "3_af2")

for directory in [rfdiff_out_dir, mpnn_out_dir, af2_out_dir]:
    os.makedirs(directory, exist_ok=True)


#____________ RFdiffusion PREP ______________

python_path_rfdiff = '/home/tsatler/anaconda3/envs/SE3nv/bin/python'
rfdiff_out_path = os.path.join(rfdiff_out_dir, "")  #Appending empty string "" results in rfdiff_out_path ending with a correct directory separator (/ or \)

# ___________ RUN RFdiffusion ________________
"""
RFdiffusion will generate new protein structures according to the contig specified
INPUT: starting pdb, contig, number of structures to design
OUTPUT: generated structure (pdb file), metadata associated with specific run for each generated structure(trb format)

TRB file contains useful information. 
In this script data from con_hal_pdb_idx/complex_con_hal_pdb_idx and 'complex_con_ref_idx0' are used in helper functions
See RFdiffusion git for details.
"""

potentials__ = 'type:olig_contacts,weight_intra:1,weight_inter:0.06'

#NOB use backslash \"{potential}\" to keep the "" correclty since the whole thing is already in ""
ckpt='/home/tsatler/RFdif/RFdiffusion/models/Base_epoch8_ckpt.pt' 
os.system(f"{python_path_rfdiff} /home/tsatler/RFdif/RFdiffusion/scripts/run_inference.py \
          inference.output_prefix={rfdiff_out_path} \
          inference.input_pdb=/{args.original_pdb_path} \
          'contigmap.contigs={args.contig}' \
          inference.num_designs={args.num_designs_rfdiff} \
          inference.symmetry='C11' \
          'potentials.guiding_potentials=[\"{potentials__}\"]' \
          potentials.olig_intra_all=True \
          potentials.olig_inter_all=True \
          potentials.guide_scale=2 \
          potentials.guide_decay='quadratic' \
          inference.ckpt_override_path={ckpt}")


"""
--config-name=symmetry \
inference.symmetry="C6" \
inference.num_designs=10 \
inference.output_prefix="example_outputs/C6_oligo" \
'potentials.guiding_potentials=["type:olig_contacts,weight_intra:1,weight_inter:0.1"]' \
potentials.olig_intra_all=True \
potentials.olig_inter_all=True \
potentials.guide_scale=2.0 \
potentials.guide_decay="quadratic" \
'contigmap.contigs=[480-480]'
"""

# ____________ MPNN PREP1 _________________
python_path_mpnn = '/home/tsatler/anaconda3/envs/mlfold/bin/python'
pymol_python_path = '/home/aljubetic/conda/envs/pyro/bin/python'

#___________MPNN PREP2_____ RECHAINING RFdiff PDBs____________
"""
RFdiffusion joins chain sequences together
For contig "[A1-30/4-6/C1-30/0 D1-30/0 B1-30]" you get two chains.
This is problematic because AF2 than folds this incorrectly as if D and B are also connected.
To solve this chain IDs are changed using rechain.py. 
The script finds chainbreaks according to pyhisical distance between CA atoms.
""" 


rf_pdbs = glob.glob(os.path.join(rfdiff_out_path, '*.pdb'))
for pdb in rf_pdbs:
    os.system(f'{pymol_python_path} /home/nbizjak/projects/11_04_2023_rigid_connections/rechain.py {pdb} {pdb}')

path_for_parsed_chains = os.path.join(mpnn_out_dir, "parsed_pdbs.jsonl")
path_for_assigned_chains = os.path.join(mpnn_out_dir, "assigned_pdbs.jsonl")
path_for_fixed_positions = os.path.join(mpnn_out_dir, "fixed_pdbs.jsonl")


#_________ MPNN PREP3______ Helper scripts _______
"""
Sripts create assigned_pdbs.jsonl, fixed_pdbs.jsonl and parsed_pdbs.jsonl which are used by ProteinMPNN
"""
os.system(f'{python_path_mpnn} /home/tsatler/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
           --input_path={rfdiff_out_dir} \
           --output_path={path_for_parsed_chains}')

os.system(f"{python_path_mpnn} /home/tsatler/ProteinMPNN/helper_scripts/assign_fixed_chains.py \
          --input_path={path_for_parsed_chains} \
          --output_path={path_for_assigned_chains} \
          --chain_list '{chains_to_design}'")


# Very important in MPNN: The first amino acid in the chain corresponds to 1 and not PDB residues index for now.
#If this changes the process_pdb_files function must change accordingly
fixed_pos_path = help_functions.process_pdb_files(rfdiff_out_dir, mpnn_out_dir) 

#_____________ RUN ProteinMPNN_____________

os.system(f'{python_path_mpnn} /home/tsatler/ProteinMPNN/protein_mpnn_run.py \
           --jsonl_path {path_for_parsed_chains} \
           --fixed_positions_jsonl {path_for_fixed_positions} \
          --chain_id_jsonl {path_for_assigned_chains} \
           --out_folder {mpnn_out_dir} \
           --num_seq_per_target {args.num_seq_per_target_mpnn} \
           --sampling_temp "0.1" \
           --batch_size 1')


#_________ AF2 PREP____ Changing MPNN fasta for AF2_________
"""
In order for AF2 to identifiy diffrent chains ":" must be placed between sequences
This is done by change_sequence_in_fasta

A specific RFdifff pdb can be provided (_0.pdb) since all RFdiff pdbs in a run have same chains
Only difference in the designed chain (chain A) which is not used by the function.
"""
# TO DO fileter duplicates
    #test ali so duplikati na nivoju mpnn
fasta_dir = os.path.join(mpnn_out_dir, "seqs")
fasta_files = sorted(glob.glob(os.path.join(fasta_dir, "*.fa"))) #glob ni sorted bo deafultu 
print(fasta_files)
rfdiff_pdb = os.path.join(rfdiff_out_path, '_0.pdb')
python_path_af2 = "/home/aljubetic/AF2/CF2.3/colabfold-conda/bin/python" # source /home/aljubetic/bin/set_up_AF2.3.sh


#________________ RUN AF2______________
for fasta_file in fasta_files: 
    model_num = help_functions.get_token_value(os.path.basename(fasta_file), "_", "(\d+)") #get 0 from _0.fa using reg exp
    model_dir = os.path.join(af2_out_dir, f"model_{model_num}") #create name for af2 directory name: model_0
    help_functions.change_sequence_in_fasta_complex(rfdiff_pdb, fasta_file)
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    os.system(f'source /home/aljubetic/bin/set_up_AF2.3.sh && {python_path_af2} /home/aljubetic/AF2/CF2.3/colabfold-conda/bin/colabfold_batch \
            --model-type alphafold2_multimer_v3 \
            --msa-mode single_sequence \
            {fasta_file} {model_dir}')

#FINAL OPERATIONS


json_directories = glob.glob(os.path.join(af2_out_dir, "*"))


for model_i in json_directories:  # for model_i in [model_0, model_1, model_2 ,...]
    
    trb_num = help_functions.get_token_value(os.path.basename(model_i), "model_", "(\d+)") #get 0 from model_0 using reg exp
    
    help_functions.rename_pdb_create_csv(args.output_dir, rfdiff_out_dir, trb_num, model_i, pdb_path)
    
    python_path = "/home/aljubetic/conda/envs/pyro/bin/python"
    
csv_path = os.path.join(args.output_dir, "output.csv") #constructed path 'output.csv defined in rename_pdb_create_csv function
os.system(f'{python_path} /home/nbizjak/projects/11_04_2023_rigid_connections/scoring_rg_charge_sap.py \
            {csv_path}')
    
scores_rg_path = os.path.join(args.output_dir, "scores_rg_charge_sap.csv") #'scores_rg_charge_sap.csv defined in scoring_rg_... script
help_functions.merge_csv(args.output_dir, csv_path, scores_rg_path)

os.remove(csv_path)
os.remove(scores_rg_path)
    
#renamed_pdb = os.path.join(os.path.dirname(args.output_dir), "final_pdbs", "*.pdb")
#df = help_functions.create_dataframe(renamed_pdb, args.output_dir)

"""
python_path = "/home/tsatler/anaconda3/envs/pyro/bin/python"
os.system(f'{python_path} /home/nbizjak/projects/11_04_2023_rigid_connections/scoring_rg_charge_sap.py \
          /home/nbizjak/projects/11_04_2023_rigid_connections/output.csv')

help_functions.merge_csv("output.csv", "scores_rg_charge_sap.csv")

os.remove("output.csv")
os.remove("scores_rg_charge_sap.csv")
"""


"""
for directory in [rfdiff_out_dir, mpnn_out_dir, af2_out_dir]:
    try:
        shutil.rmtree(directory)
        print(f"Directory '{directory}' has been deleted successfully!")
    except OSError as e:
        print(f"Error deleting the directory '{directory}': {e}")
"""
