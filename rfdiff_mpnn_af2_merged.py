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

def run_and_log(command, log_func=log.info, dry_run=False):
    """Runs a command using os.system and also logs the command before running using log.info"""
    if log_func:
        log_func(command)
    if not dry_run:
        os.system(command)

scripts_folder = pathlib.Path(__file__).resolve().parent / "scripts"



def general_config_prep(cfg):
    log.info("Running generalPrep")
    os.makedirs(cfg.output_dir, exist_ok=True)

    # We want to have all the params in the cfg struct. Thus, use we open_dict to be able to write new data.
    with open_dict(cfg):
        cfg.output_dir = str(pathlib.Path(cfg.output_dir).resolve()) # Output dir MUST be absolute, otherwise ProtMPPN complains about missing seq_chain (update: it probably doesn't matter)
        cfg.rfdiff_out_dir = os.path.join(cfg.output_dir, "1_rfdiff")
        cfg.mpnn_out_dir = os.path.join(cfg.output_dir, "2_mpnn")
        cfg.af2_out_dir = os.path.join(cfg.output_dir, "3_af2")
        cfg.rfdiff_out_path = os.path.join(cfg.rfdiff_out_dir, "")  #Appending empty string "" results in rfdiff_out_path ending with a correct directory separator (/ or \)

        cfg.path_for_parsed_chains = os.path.join(cfg.mpnn_out_dir, "parsed_pdbs.jsonl")
        cfg.path_for_assigned_chains = os.path.join(cfg.mpnn_out_dir, "assigned_pdbs.jsonl")
        cfg.path_for_fixed_positions = os.path.join(cfg.mpnn_out_dir, "fixed_pdbs.jsonl")

        cfg.fasta_dir = os.path.join(cfg.mpnn_out_dir, "seqs")
        cfg.rfdiff_pdb = os.path.join(cfg.rfdiff_out_path, '_0.pdb')

        if cfg.get("skipRfDiff", False):
            # We only need to redesign the chains specified in designable_residues
            cfg.chains_to_design = " ".join(sorted({_[0] for _ in cfg.designable_residues}))
            log.info(f"Skipping RFdiff, only redesigning chains specified in designable_residues: {cfg.chains_to_design}")

        # I suggest the following: count("/0", contig) -> chains_to_design = " ".join(chain_letters[:count]), unless specified (in run.yaml, it should be null, None or sth similar)
        chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" # What happens after 26 chains? RfDiff only supports up to 26 chains: https://github.com/RosettaCommons/RFdiffusion/blob/ba8446eae0fb80c121829a67d3464772cc827f01/rfdiffusion/contigs.py#L40C29-L40C55
        if cfg.chains_to_design == None: #TODO this will likely break for sym mode
            breaks = cfg.contig.count("/0 ") + 1
            cfg.chains_to_design = ' '.join(chain_letters[:breaks])
            log.info(f"Chains to design (according to contig chain breaks): {cfg.chains_to_design}")
            



    
    for directory in [cfg.rfdiff_out_dir, cfg.mpnn_out_dir, cfg.af2_out_dir]:
        os.makedirs(directory, exist_ok=True)
    
    # No need to return; cfg is mutable.
    

def run_rfdiff(cfg):
    """
    RFdiffusion will generate new protein structures according to the contig specified
    INPUT: starting pdb, contig, number of structures to design
    OUTPUT: generated structure (pdb file), metadata associated with specific run for each generated structure (trb format)

    TRB file contains useful information. 
    In this script data from con_hal_pdb_idx/complex_con_hal_pdb_idx and 'complex_con_ref_idx0' are used in helper functions
    See RFdiffusion git for details.
    """
    log.info("***************Running runRFdiff***************")
    rfdiff_cmd_str = f"{cfg.python_path_rfdiff} {cfg.inference_path_rfdiff} \
          inference.output_prefix={cfg.rfdiff_out_path} \
          'contigmap.contigs={cfg.contig}' \
          inference.num_designs={cfg.num_designs_rfdiff} \
          -cn prosculpt2rfdiff.yaml -cd {cfg.output_dir}"
    run_and_log(rfdiff_cmd_str)
    
    log.info("***************After running RFdiffusion***************")

def rechain_rfdiff_pdbs(cfg):
    """
    RFdiffusion joins chain sequences together
    For contig "[A1-30/4-6/C1-30/0 D1-30/0 B1-30]" you get two chains.
    This is problematic because AF2 than folds this incorrectly as though D and B were also connected.
    To solve this chain IDs are changed using rechain.py. 
    The script finds chainbreaks according to physical distance between CA atoms.
    """ 
    log.info("Running rechainRFdiffPDBs")
    rf_pdbs = glob.glob(os.path.join(cfg.rfdiff_out_path, '*.pdb'))
    for pdb in rf_pdbs:
        run_and_log(
            f'{cfg.pymol_python_path} {scripts_folder / "rechain.py"} {pdb} {pdb} --chain_break_cutoff_A {cfg.chain_break_cutoff_A}'
        )
        

    log.info("After rechaining")

def parse_additional_args(cfg, group):
    dodatniArgumenti = ""
    for k, v in (cfg.get(group, {}) or {}).items(): # or to allow for empty groups
        dodatniArgumenti += f" {k} {v}"
    return(dodatniArgumenti)

def do_cycling(cfg):
    log.info("Running do_cycling")
    for cycle in range(cfg.af2_mpnn_cycles):
        print("cycleeeeee", cycle)

        trb_paths = None
        input_mpnn = cfg.rfdiff_out_dir # First cycle has this, other cycles overwrite this var

        if not cycle == 0: # All other cycles get starting PDBs from AF2
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
        
        run_and_log(
            f'{cfg.python_path_mpnn} {os.path.join(cfg.mpnn_installation_path, "helper_scripts", "parse_multiple_chains.py")} \
            --input_path={input_mpnn} \
            --output_path={cfg.path_for_parsed_chains}'       
        )


        run_and_log(
            f"{cfg.python_path_mpnn} {os.path.join(cfg.mpnn_installation_path, 'helper_scripts', 'assign_fixed_chains.py')} \
                --input_path={cfg.path_for_parsed_chains} \
                --output_path={cfg.path_for_assigned_chains} \
                --chain_list='{cfg.chains_to_design}'"
                )

        fixed_pos_path = prosculpt.process_pdb_files(input_mpnn, cfg.mpnn_out_dir, cfg, trb_paths) # trb_paths is optional (default: None) and only used in non-first cycles
        # trb_paths is atm not used in process_pdb_files anyway -- a different approach is used (file.pdb -> withExtension .trb), which ensures the PDB and TRB files match.
        log.info(f"Fixed positions path: {fixed_pos_path}")

        #_____________ RUN ProteinMPNN_____________
        # At first cycle, use num_seq_per_target from config. In subsequent cycles, set it to 1.
        proteinMPNN_cmd_str = f'{cfg.python_path_mpnn} {os.path.join(cfg.mpnn_installation_path, "protein_mpnn_run.py")} \
            --jsonl_path {cfg.path_for_parsed_chains} \
            --fixed_positions_jsonl {cfg.path_for_fixed_positions} \
            --chain_id_jsonl {cfg.path_for_assigned_chains} \
            --out_folder {cfg.mpnn_out_dir} \
            --num_seq_per_target {cfg.num_seq_per_target_mpnn if cycle == 0 else 1} \
            {"--sampling_temp 0.1 --backbone_noise 0" if cfg.get("skipRfDiff", False) else ""} \
            {parse_additional_args(cfg, "pass_to_mpnn")} \
            --batch_size 1'
        
        run_and_log(proteinMPNN_cmd_str)

        log.info("Preparing to empty af2 directory.")
        
        # af2 directory must be empty
        shutil.rmtree(cfg.af2_out_dir)
        os.makedirs(cfg.af2_out_dir, exist_ok=True)


        #________________ RUN AF2______________
        fasta_files = sorted(glob.glob(os.path.join(cfg.fasta_dir, "*.fa"))) # glob is not sorted by default
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

            model_order = str(cfg.model_order).split(",") # If only one number is passed, Hydra converts it to int
            num_models = len(model_order) #model order "1,2,3" -> num of models = 3
            if cycle == 0: #have to run af2 differently in first cycle 
                run_and_log(
                    f'source {cfg.af_setup_path} && {cfg.python_path_af2} {cfg.colabfold_setup_path} \
                                --model-type alphafold2_multimer_v3 \
                                --msa-mode single_sequence \
                                {fasta_file} {model_dir} \
                                --model-order {cfg.model_order} \
                                {parse_additional_args(cfg, "pass_to_af")} \
                                --num-models {num_models}'
                                )
            else: 
                for model_number in model_order:
                    if af2_model_num == model_number: # From the af2 model 4 you want only model 4 not also 2 and for 2 only 2 not 4 (--model_order "2,4")
                        num_models = 1
                        run_and_log(
                            f'source {cfg.af_setup_path} && {cfg.python_path_af2} {cfg.colabfold_setup_path} \
                                --model-type alphafold2_multimer_v3 \
                                --msa-mode single_sequence \
                                {fasta_file} {model_dir} \
                                --model-order {model_number} \
                                {parse_additional_args(cfg, "pass_to_af")} \
                                --num-models {num_models}'
                                )

        # msa single sequence makes sense for designed proteins

def final_operations(cfg):
    log.info("Final operations")
    json_directories = glob.glob(os.path.join(cfg.af2_out_dir, "*"))

    for model_i in json_directories:  # for model_i in [model_0, model_1, model_2 ,...]
        
        trb_num = prosculpt.get_token_value(os.path.basename(model_i), "model_", "(\d+)") #get 0 from model_0 using reg exp
        prosculpt.rename_pdb_create_csv(cfg.output_dir, cfg.rfdiff_out_dir, trb_num, model_i, cfg.pdb_path)
        
                
    csv_path = os.path.join(cfg.output_dir, "output.csv") #constructed path 'output.csv defined in rename_pdb_create_csv function
    run_and_log(
        f'{cfg.python_path} {scripts_folder / "scoring_rg_charge_sap.py"} {csv_path}'
        )
        
    scores_rg_path = os.path.join(cfg.output_dir, "scores_rg_charge_sap.csv") #'scores_rg_charge_sap.csv defined in scoring_rg_... script
    prosculpt.merge_csv(cfg.output_dir, csv_path, scores_rg_path)

    os.remove(csv_path)
    os.remove(scores_rg_path)


def pass_config_to_rfdiff(cfg):
    """
    This function creates a new .yaml file with "inference", "potentials" etc. read from run.yaml. Also add defaults: - base \n - _self_
    Then, when calling RfDiff, add -cd [folder_with_this_created_config] -cn [name_of_created_config_file]
    """
    # Keep in mind that if passing through cmd, these parameters need to have ++ in front of them.
    # TODO: should everything go through the command line? Makes it a bit cleaner.
    def keep_only_groups(pair):
        groupsToPassToRfDiff = cfg.pass_to_rfdiff # groups that should be passed to RfDiff
        key, value = pair
        if key not in groupsToPassToRfDiff:
            return False
        else:
            return True


    new_cfg = dict(filter(keep_only_groups, cfg.items()))
    new_cfg["defaults"] = ["base", "_self_"] # TODO: is the base always the right thing to include? For symm it should be symmetry.
    with open_dict(new_cfg["inference"]):
        new_cfg["inference"]["input_pdb"] = cfg.pdb_path

    log.info(f"Saving this config for RfDiff: {OmegaConf.to_yaml(new_cfg)}")
    # save to file
    with open(os.path.join(cfg.output_dir, "prosculpt2rfdiff.yaml"), "w") as f:
        f.write(OmegaConf.to_yaml(new_cfg))


    pass

@hydra.main(version_base=None, config_path="config", config_name="run")
def prosculptApp(cfg: DictConfig) -> None:
    log.info("Hydrić")
    log.info("The following configuration was passed: \n" + OmegaConf.to_yaml(cfg))

    #log.debug(f"Now in Hydra, cwd = {os.getcwd()}")

    general_config_prep(cfg)

    if not cfg.get("skipRfDiff", False):
        pass_config_to_rfdiff(cfg)
        run_rfdiff(cfg)
        rechain_rfdiff_pdbs(cfg)
    else:
        log.info("*** Skipping RfDiff ***")
        # Copy input PDB to RfDiff_output_dir and rename it to follow the token scheme
        shutil.copy(cfg.pdb_path, os.path.join(cfg.rfdiff_out_dir, "_0.pdb"))

    do_cycling(cfg)
    final_operations(cfg)


if __name__ == "__main__":
    print("File: ", __file__)
    a = pathlib.Path(__file__).resolve().parent
    log.info(f"Before running MYapp, cwd = {os.getcwd()}")
    prosculptApp()