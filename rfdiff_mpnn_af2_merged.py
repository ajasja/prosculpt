from omegaconf import DictConfig, OmegaConf
import hydra
import os
from hydra.utils import get_original_cwd, to_absolute_path
import logging
import prosculpt
import glob
import re
import shutil
import pathlib
from omegaconf import open_dict




log = logging.getLogger(__name__)
scripts_folder = pathlib.Path(__file__).resolve().parent / "scripts"

def generalPrep(cfg):
    log.info("Running generalPrep")
    os.makedirs(cfg.output_dir, exist_ok=True)

    # We want to have all the params in the cfg struct. Thus, use we open_dict to be able to write new data.
    with open_dict(cfg):
        cfg.output_dir = str(pathlib.Path(cfg.output_dir).resolve()) # Output dir MUST be absolute, otherwise ProtMPPN complains about missing seq_chain (upadte: it probably doesn't matter)
        cfg.rfdiff_out_dir = os.path.join(cfg.output_dir, "1_rfdiff")
        cfg.mpnn_out_dir = os.path.join(cfg.output_dir, "2_mpnn")
        cfg.af2_out_dir = os.path.join(cfg.output_dir, "3_af2")
        cfg.rfdiff_out_path = os.path.join(cfg.rfdiff_out_dir, "")  #Appending empty string "" results in rfdiff_out_path ending with a correct directory separator (/ or \)

        cfg.path_for_parsed_chains = os.path.join(cfg.mpnn_out_dir, "parsed_pdbs.jsonl")
        cfg.path_for_assigned_chains = os.path.join(cfg.mpnn_out_dir, "assigned_pdbs.jsonl")
        cfg.path_for_fixed_positions = os.path.join(cfg.mpnn_out_dir, "fixed_pdbs.jsonl")

        cfg.fasta_dir = os.path.join(cfg.mpnn_out_dir, "seqs")
        cfg.rfdiff_pdb = os.path.join(cfg.rfdiff_out_path, '_0.pdb')
        python_path_af2 = "/home/aljubetic/AF2/CF2.3/colabfold-conda/bin/python" # source /home/aljubetic/bin/set_up_AF2.3.sh


    
    for directory in [cfg.rfdiff_out_dir, cfg.mpnn_out_dir, cfg.af2_out_dir]:
        os.makedirs(directory, exist_ok=True)
    
    log.info(cfg)
    # No need to return; cfg is mutable.
    #return(rfdiff_out_dir, mpnn_out_dir, af2_out_dir, rfdiff_out_path)

def runRFdiff(cfg):
    """
    RFdiffusion will generate new protein structures according to the contig specified
    INPUT: starting pdb, contig, number of structures to design
    OUTPUT: generated structure (pdb file), metadata associated with specific run for each generated structure (trb format)

    TRB file contains useful information. 
    In this script data from con_hal_pdb_idx/complex_con_hal_pdb_idx and 'complex_con_ref_idx0' are used in helper functions
    See RFdiffusion git for details.
    """
    log.info("Running runRFdiff")

    log.info(f"{cfg.python_path_rfdiff} {cfg.inference_path_rfdiff} \
          inference.output_prefix={cfg.rfdiff_out_path} \
          inference.input_pdb={cfg.pdb_path} \
          'contigmap.contigs={cfg.contig}' \
          inference.num_designs={cfg.num_designs_rfdiff}")
    
    os.system(f"{cfg.python_path_rfdiff} {cfg.inference_path_rfdiff} \
          inference.output_prefix={cfg.rfdiff_out_path} \
          inference.input_pdb={cfg.pdb_path} \
          'contigmap.contigs={cfg.contig}' \
          inference.num_designs={cfg.num_designs_rfdiff}")
    
    log.info("After running RFdiffusion")

def rechainRFdiffPDBs(cfg):
    """
    RFdiffusion joins chain sequences together
    For contig "[A1-30/4-6/C1-30/0 D1-30/0 B1-30]" you get two chains.
    This is problematic because AF2 than folds this incorrectly as though D and B were also connected.
    To solve this chain IDs are changed using rechain.py. 
    The script finds chainbreaks according to pyhisical distance between CA atoms.
    """ 
    log.info("Running rechainRFdiffPDBs")
    rf_pdbs = glob.glob(os.path.join(cfg.rfdiff_out_path, '*.pdb'))
    for pdb in rf_pdbs:

        log.info(f'{cfg.pymol_python_path} {scripts_folder / "rechain.py"} {pdb} {pdb} --chain_break_cutoff_A {cfg.chain_break_cutoff_A}')

        os.system(f'{cfg.pymol_python_path} {scripts_folder / "rechain.py"} {pdb} {pdb} --chain_break_cutoff_A {cfg.chain_break_cutoff_A}')
    log.info("After rechaining")

def do_cycling(cfg):
    log.info("Running do_cycling")
    log.info("Entering the spinning loop of cycles")
    for cycle in range(cfg.af2_mpnn_cycles):
        print("cycleeeeee", cycle)

        trb_paths = None
        input_mpnn = cfg.rfdiff_out_dir # First cycle has this, other cycles overwrite this var

        if not cycle == 0: # All other cycles get starting PDBs from AF2
            print("______________________________entering cycle_____________________________________")
            cycle_directory = os.path.join(cfg.output_dir, "2_1_cycle_directory")
            if os.path.exists(cycle_directory):
                shutil.rmtree(cycle_directory)
            os.makedirs(cycle_directory, exist_ok=True)

            af2_model_subdicts = glob.glob(os.path.join(cfg.af2_out_dir, "*"))

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

            shutil.rmtree(cfg.mpnn_out_dir) # Remove MPNN dir so you can create new sequences
            os.makedirs(cfg.mpnn_out_dir, exist_ok=True) # Create it again
            trb_paths = os.path.join(cfg.rfdiff_out_dir, "*.trb")
        #endif
        # All cycles run the same commands
        log.info(f'{cfg.python_path_mpnn} {os.path.join(cfg.mpnn_installation_path, "helper_scripts", "parse_multiple_chains.py")} \
            --input_path={input_mpnn} \
            --output_path={cfg.path_for_parsed_chains}')

        os.system(f'{cfg.python_path_mpnn} {os.path.join(cfg.mpnn_installation_path, "helper_scripts", "parse_multiple_chains.py")} \
            --input_path={input_mpnn} \
            --output_path={cfg.path_for_parsed_chains}')

        log.info(f"{cfg.python_path_mpnn} {os.path.join(cfg.mpnn_installation_path, 'helper_scripts', 'assign_fixed_chains.py')} \
                --input_path={cfg.path_for_parsed_chains} \
                --output_path={cfg.path_for_assigned_chains} \
                --chain_list='{cfg.chains_to_design}'")

        os.system(f"{cfg.python_path_mpnn} {os.path.join(cfg.mpnn_installation_path, 'helper_scripts', 'assign_fixed_chains.py')} \
                --input_path={cfg.path_for_parsed_chains} \
                --output_path={cfg.path_for_assigned_chains} \
                --chain_list='{cfg.chains_to_design}'")

        fixed_pos_path = prosculpt.process_pdb_files(input_mpnn, cfg.mpnn_out_dir, trb_paths) # trb_paths is optional (default: None) and only used in non-first cycles
        log.info(f"Fixed positions path: {fixed_pos_path}")

        #_____________ RUN ProteinMPNN_____________
        # At first cycle, use num_seq_per_target from config. In subsequent cycles, set it to 1.
        log.info(f'{cfg.python_path_mpnn} /home/tsatler/ProteinMPNN/protein_mpnn_run.py \
            --jsonl_path {cfg.path_for_parsed_chains} \
            --fixed_positions_jsonl {cfg.path_for_fixed_positions} \
            --chain_id_jsonl {cfg.path_for_assigned_chains} \
            --out_folder {cfg.mpnn_out_dir} \
            --num_seq_per_target {cfg.num_seq_per_target_mpnn if cycle == 0 else 1} \
            --sampling_temp "0.1" \
            --batch_size 1')

        os.system(f'{cfg.python_path_mpnn} /home/tsatler/ProteinMPNN/protein_mpnn_run.py \
            --jsonl_path {cfg.path_for_parsed_chains} \
            --fixed_positions_jsonl {cfg.path_for_fixed_positions} \
            --chain_id_jsonl {cfg.path_for_assigned_chains} \
            --out_folder {cfg.mpnn_out_dir} \
            --num_seq_per_target {cfg.num_seq_per_target_mpnn if cycle == 0 else 1} \
            --sampling_temp "0.1" \
            --batch_size 1')

        log.info("Preparing to empty af2 directory.")
        
        # af2 directory must be empty
        shutil.rmtree(cfg.af2_out_dir)
        os.makedirs(cfg.af2_out_dir, exist_ok=True)


        #________________ RUN AF2______________
        fasta_files = sorted(glob.glob(os.path.join(cfg.fasta_dir, "*.fa"))) #glob ni sorted po deafultu 
        for fasta_file in fasta_files: # Number of fasta files corresponds to number of rfdiff models
            if cycle == 0:
                rf_model_num = prosculpt.get_token_value(os.path.basename(fasta_file), "_", "(\d+)") #get 0 from _0.fa using reg exp
            else:
                rf_model_num = prosculpt.get_token_value(os.path.basename(fasta_file), "rf_", "(\d+)") #get 0 from rf_0__model_1__cycle_2__itr_0__.pdb
                af2_model_num = prosculpt.get_token_value(os.path.basename(fasta_file), "model_", "(\d+)") #get 1 from rf_0__model_1__cycle_2__itr_0__.pdb
            model_dir = os.path.join(cfg.af2_out_dir, f"model_{rf_model_num}") #create name for af2 directory name: model_0
            prosculpt.change_sequence_in_fasta(cfg.rfdiff_pdb, fasta_file)
            if not os.path.exists(model_dir):
                os.makedirs(model_dir)

            model_order = cfg.model_order.split(",")
            num_models = len(model_order) #model order "1,2,3" -> num of models = 3
            if cycle == 0: #have to run af2 differently in first cycle 
                os.system(f'source {cfg.af_setup_path} && {cfg.python_path_af2} {cfg.colabfold_setup_path} \
                                --model-type alphafold2_multimer_v3 \
                                --msa-mode single_sequence \
                                {fasta_file} {model_dir} \
                                --model-order {cfg.model_order} \
                                --num-models {num_models}')
            else: 
                for model_number in model_order:
                    if af2_model_num == model_number: # From the af2 model 4 you want only model 4 not also 2 and for 2 only 2 not 4 (--model_order "2,4")
                        num_models = 1
                        os.system(f'source {cfg.af_setup_path} && {cfg.python_path_af2} {cfg.colabfold_setup_path} \
                                --model-type alphafold2_multimer_v3 \
                                --msa-mode single_sequence \
                                {fasta_file} {model_dir} \
                                --model-order {model_number} \
                                --num-models {num_models}')

        # msa single sequence makes sense for designed (no sense to aligbn to natural proteins)

def finalOperations(cfg):
    log.info("Final operations")
    json_directories = glob.glob(os.path.join(cfg.af2_out_dir, "*"))

    for model_i in json_directories:  # for model_i in [model_0, model_1, model_2 ,...]
        
        trb_num = prosculpt.get_token_value(os.path.basename(model_i), "model_", "(\d+)") #get 0 from model_0 using reg exp
        prosculpt.rename_pdb_create_csv(cfg.output_dir, cfg.rfdiff_out_dir, trb_num, model_i, cfg.pdb_path)
        
                
    csv_path = os.path.join(cfg.output_dir, "output.csv") #constructed path 'output.csv defined in rename_pdb_create_csv function
    os.system(f'{cfg.python_path} {scripts_folder / "scoring_rg_charge_sap.py"} \
                {csv_path}')
        
    scores_rg_path = os.path.join(cfg.output_dir, "scores_rg_charge_sap.csv") #'scores_rg_charge_sap.csv defined in scoring_rg_... script
    prosculpt.merge_csv(cfg.output_dir, csv_path, scores_rg_path)

    os.remove(csv_path)
    os.remove(scores_rg_path)


@hydra.main(version_base=None, config_path="config", config_name="run")
def meinApp(cfg: DictConfig) -> None:
    log.info("The following config was passed: ")
    log.info(cfg)
    log.info("Hydrić")
    log.info(OmegaConf.to_yaml(cfg))

    log.info(f"Now in Hydra, cwd = {os.getcwd()}")

    if cfg.chains_to_design != "A":
        log.warn("chains_to_design is set to something else than A. Please note that RfDiff sometimes renames all chains to A. Please check the log for potential errors.")
        log.info(f"Chains to design: {cfg.chains_to_design}")
    else:
        log.info("Chains to design: A")
        log.info(cfg.chains_to_design)
    
    generalPrep(cfg)
    runRFdiff(cfg)
    rechainRFdiffPDBs(cfg)
    do_cycling(cfg)
    finalOperations(cfg)

#TODO: Actually merge with _simmetry

if __name__ == "__main__":
    print("File: ", __file__)
    a = pathlib.Path(__file__).resolve().parent
    print(a)
    print(f"Before runnning MYapp, cwd = {os.getcwd()}")
    meinApp()