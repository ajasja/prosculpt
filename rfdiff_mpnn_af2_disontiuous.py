import os
import prosculpt
import glob
import argparse
import re
import shutil



#/home/nbizjak/projects/11_04_2023_rigid_connections/inputs_for_rfdiff/phosphoCC_bundle.pdb
#/home/nbizjak/projects/11_04_2023_rigid_connections/inputs_for_rfdiff/dva_helixa_delujoca.pdb
pdb_path ="/home/nbizjak/projects/11_04_2023_rigid_connections/inputs_for_rfdiff/mALb8-antiA.pdb"
output_dir = "/home/nbizjak/prosculpt_dev/outputs/deafult_output/"
contig = "[C33-60/5-7/A1-30/0 B61-120]" #[A1-30/5/C1-30/0 D1-30/0  B1-30] #[A1-37/30-60/A42-70] dva helixa  #[C33-60/3/A1-30/0 B61-120] mlba... [C33-60/5/A1-30/0 B61-120]
num_designs_rfdiff = 1
num_seq_per_target_mpnn = 1
chains_to_design = "A" #might be very important for complexes


#if designing complexes change in helper funcitions in the pdb proceisng function

parser = argparse.ArgumentParser(description="Run protein design pipeline")
parser.add_argument("--contig", default=contig, help="Contig to design")
parser.add_argument("--original_pdb_path", default=pdb_path, help="Path to input PDB file")
parser.add_argument("--output_dir", default = output_dir, help="Path to output directory")
parser.add_argument("--num_designs_rfdiff", type=int, default=num_designs_rfdiff, help="Number of designs to generate with RFdiffusion")
parser.add_argument("--num_seq_per_target_mpnn", type=int, default=num_seq_per_target_mpnn, help="Number of sequences to generate per target with MPNN")
parser.add_argument("--chains_to_design_mpnn", default= chains_to_design, help="All chains of protein even if not redesigned (except if doing binder design)") 
parser.add_argument("--af2_mpnn_cycles", default = 1, type = int, help="Cycles of mpnn-af2. Deafult is 1")
parser.add_argument("--model_order", default = "1,2,3,4,5",type = str, help="Choose specific models to cycle in af2-mpnn")


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

os.system(f"{python_path_rfdiff} /home/tsatler/RFdif/RFdiffusion/scripts/run_inference.py \
          inference.output_prefix={rfdiff_out_path} \
          inference.input_pdb=/{args.original_pdb_path} \
          'contigmap.contigs={args.contig}' \
          inference.num_designs={args.num_designs_rfdiff}")



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
    os.system(f'{pymol_python_path} /home/nbizjak/prosculpt_dev/prosculpt/rechain.py {pdb} {pdb}')

print("entering for loop of cycles")
#CYCLES OF MPNN AND AF2
for cycle in range(args.af2_mpnn_cycles):

    print("cycleeeeee", cycle)
    path_for_parsed_chains = os.path.join(mpnn_out_dir, "parsed_pdbs.jsonl")
    path_for_assigned_chains = os.path.join(mpnn_out_dir, "assigned_pdbs.jsonl")
    path_for_fixed_positions = os.path.join(mpnn_out_dir, "fixed_pdbs.jsonl")
    
    if cycle == 0: #For the first cycle...
        input_mpnn = rfdiff_out_dir

        os.system(f'{python_path_mpnn} /home/tsatler/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
            --input_path={input_mpnn} \
            --output_path={path_for_parsed_chains}')

        os.system(f"{python_path_mpnn} /home/tsatler/ProteinMPNN/helper_scripts/assign_fixed_chains.py \
                --input_path={path_for_parsed_chains} \
                --output_path={path_for_assigned_chains} \
                --chain_list '{chains_to_design}'")

        fixed_pos_path = prosculpt.process_pdb_files(input_mpnn, mpnn_out_dir)
    else: #All other cycles get starting pdbs from af2
        print("______________________________entering cycle_____________________________________")
        cycle_directory = os.path.join(args.output_dir, "2_1_cycle_directory")
        if os.path.exists(cycle_directory):
           shutil.rmtree(cycle_directory)
        os.makedirs(cycle_directory, exist_ok=True)

        af2_model_subdicts = glob.glob(os.path.join(af2_out_dir, "*"))

        #Access all af2 models and put them in one intemediate directory to get the same ored as in the 1st cycle (all pdbs in one directory)
        for model_subdict in af2_model_subdicts: 
            af2_pdbs = sorted(glob.glob(os.path.join(model_subdict, "T*.pdb")))
            for i, af2_pdb in enumerate(af2_pdbs):
                
                af_model_num = prosculpt.get_token_value(os.path.basename(af2_pdb), "model_", "(\d+)")
                #if str(af_model_num) in args.af2_models:
                # Rename pdbs to keep the model_num traceability with orginal rfdiff structure and enable filtering which models for next cycle
                if cycle == 1:
                    rf_model_num = prosculpt.get_token_value(os.path.basename(model_subdict), "model_", "(\d+)")
                    #else:
                        #rf_model_num = prosculpt.get_token_value(os.path.basename(af2_pdb), "rf__", "(\d+)") #after first cycling modelXX directories in af_output do not correspond to rf model anymore
                shutil.move(af2_pdb, os.path.join(cycle_directory, f"rf_{rf_model_num}__model_{af_model_num}__cycle_{cycle}__itr_{i}__.pdb")) 
                            #rf_ --> rfdiffusion structure number (in rfdiff outou dir)
                            #model_ -> af2 model num, used for filtering which to cycle (preference for model 4)
                            #itr_ -> to differentiate models in later cycles (5 pdbs for model 4 from rf 0 for example)
                            # is it maybe possible to filter best ranked by af2 from the itr numbers?

        input_mpnn = cycle_directory

        shutil.rmtree(mpnn_out_dir) # Remove MPNN dir so you can create new sequences
        os.makedirs(mpnn_out_dir, exist_ok=True) # Create it again
        trb_path = os.path.join(rfdiff_out_dir, "_0.trb")

        os.system(f'{python_path_mpnn} /home/tsatler/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
            --input_path={input_mpnn} \
            --output_path={path_for_parsed_chains}')

        os.system(f"{python_path_mpnn} /home/tsatler/ProteinMPNN/helper_scripts/assign_fixed_chains.py \
                --input_path={path_for_parsed_chains} \
                --output_path={path_for_assigned_chains} \
                --chain_list '{chains_to_design}'")

        fixed_pos_path = prosculpt.process_pdb_files(input_mpnn, mpnn_out_dir, trb_path)

     
    #_____________ RUN ProteinMPNN_____________

    os.system(f'{python_path_mpnn} /home/tsatler/ProteinMPNN/protein_mpnn_run.py \
            --jsonl_path {path_for_parsed_chains} \
            --fixed_positions_jsonl {path_for_fixed_positions} \
            --chain_id_jsonl {path_for_assigned_chains} \
            --out_folder {mpnn_out_dir} \
            --num_seq_per_target {args.num_seq_per_target_mpnn} \
            --sampling_temp "0.1" \
            --batch_size 1')

    # af2 directory must be empty
    shutil.rmtree(af2_out_dir)
    os.makedirs(af2_out_dir, exist_ok=True)

    #_________ AF2 PREP____ Changing MPNN fasta for AF2_________

    fasta_dir = os.path.join(mpnn_out_dir, "seqs")
    fasta_files = sorted(glob.glob(os.path.join(fasta_dir, "*.fa"))) #glob ni sorted bo deafultu 
    rfdiff_pdb = os.path.join(rfdiff_out_path, '_0.pdb')
    python_path_af2 = "/home/aljubetic/AF2/CF2.3/colabfold-conda/bin/python" # source /home/aljubetic/bin/set_up_AF2.3.sh


    #________________ RUN AF2______________
    for fasta_file in fasta_files: # Number of fasta files corresponds to number of rfdiff models
        if cycle == 0:
            rf_model_num = prosculpt.get_token_value(os.path.basename(fasta_file), "_", "(\d+)") #get 0 from _0.fa using reg exp
        else:
            rf_model_num = prosculpt.get_token_value(os.path.basename(fasta_file), "rf_", "(\d+)") #get 0 from rf_0__model_1__cycle_2__itr_0__.pdb
            af2_model_num = prosculpt.get_token_value(os.path.basename(fasta_file), "model_", "(\d+)") #get 1 from rf_0__model_1__cycle_2__itr_0__.pdb
        model_dir = os.path.join(af2_out_dir, f"model_{rf_model_num}") #create name for af2 directory name: model_0
        prosculpt.change_sequence_in_fasta(rfdiff_pdb, fasta_file)
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)

        model_order = args.model_order.split(",")
        num_models = len(model_order) #model order "1,2,3" -> num of models = 3
        if cycle == 0: #have to run af2 differently in first cycle 
            os.system(f'source /home/aljubetic/bin/set_up_AF2.3.sh && {python_path_af2} /home/aljubetic/AF2/CF2.3/colabfold-conda/bin/colabfold_batch \
                            --model-type alphafold2_multimer_v3 \
                            --msa-mode single_sequence \
                            {fasta_file} {model_dir} \
                            --model-order {args.model_order} \
                            --num-models {num_models}')
        else: 
            for model_number in model_order:
                if af2_model_num == model_number: # From the af2 model 4 you want only model 4 not also 2 and for 2 only 2 not 4 (--model_order "2,4")
                    num_models = 1
                    os.system(f'source /home/aljubetic/bin/set_up_AF2.3.sh && {python_path_af2} /home/aljubetic/AF2/CF2.3/colabfold-conda/bin/colabfold_batch \
                            --model-type alphafold2_multimer_v3 \
                            --msa-mode single_sequence \
                            {fasta_file} {model_dir} \
                            --model-order {model_number} \
                            --num-models {num_models}')

    # msa single sequence makes sense for designed (no sense to aligbn to natural proteins)

#FINAL OPERATIONS


json_directories = glob.glob(os.path.join(af2_out_dir, "*"))


for model_i in json_directories:  # for model_i in [model_0, model_1, model_2 ,...]
    
    trb_num = prosculpt.get_token_value(os.path.basename(model_i), "model_", "(\d+)") #get 0 from model_0 using reg exp
    prosculpt.rename_pdb_create_csv(args.output_dir, rfdiff_out_dir, trb_num, model_i, pdb_path)
    
    
python_path = "/home/aljubetic/conda/envs/pyro/bin/python"
    
csv_path = os.path.join(args.output_dir, "output.csv") #constructed path 'output.csv defined in rename_pdb_create_csv function
os.system(f'{python_path} /home/nbizjak/prosculpt_dev/prosculpt/scoring_rg_charge_sap.py \
            {csv_path}')
    
scores_rg_path = os.path.join(args.output_dir, "scores_rg_charge_sap.csv") #'scores_rg_charge_sap.csv defined in scoring_rg_... script
prosculpt.merge_csv(args.output_dir, csv_path, scores_rg_path)

os.remove(csv_path)
os.remove(scores_rg_path)
    
#renamed_pdb = os.path.join(os.path.dirname(args.output_dir), "final_pdbs", "*.pdb")
#df = prosculpt.create_dataframe(renamed_pdb, args.output_dir)

"""
python_path = "/home/tsatler/anaconda3/envs/pyro/bin/python"
os.system(f'{python_path} /home/nbizjak/projects/11_04_2023_rigid_connections/scoring_rg_charge_sap.py \
          /home/nbizjak/projects/11_04_2023_rigid_connections/output.csv')

prosculpt.merge_csv("output.csv", "scores_rg_charge_sap.csv")

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

