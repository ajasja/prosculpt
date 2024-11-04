from omegaconf import DictConfig, OmegaConf
import hydra
from hydra.core.hydra_config import HydraConfig
import os
from hydra.utils import get_original_cwd, to_absolute_path
import logging
import prosculpt
import glob
import re
import shutil
import pathlib
from omegaconf import open_dict
from Bio import SeqIO


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
        cfg.path_for_tied_positions = os.path.join(cfg.mpnn_out_dir, "tied_pdbs.jsonl")

        cfg.fasta_dir = os.path.join(cfg.mpnn_out_dir, "seqs")
        cfg.rfdiff_pdb = os.path.join(cfg.rfdiff_out_path, '_0.pdb')

        if cfg.get("skipRfDiff", False):
            # We only need to redesign the chains specified in designable_residues
            cfg.chains_to_design = " ".join(sorted({_[0] for _ in cfg.designable_residues}))
            log.info(f"Skipping RFdiff, only redesigning chains specified in designable_residues: {cfg.chains_to_design}")

        if 'symmetry' not in cfg.inference:
            cfg.inference.symmetry=None

        if 'omit_AAs' not in cfg:
            cfg.omit_AAs="X" #This is the default also in proteinMPNN

        # I suggest the following: count("/0", contig) -> chains_to_design = " ".join(chain_letters[:count]), unless specified (in run.yaml, it should be null, None or sth similar)
        chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" # What happens after 26 chains? RfDiff only supports up to 26 chains: https://github.com/RosettaCommons/RFdiffusion/blob/ba8446eae0fb80c121829a67d3464772cc827f01/rfdiffusion/contigs.py#L40C29-L40C55
        if cfg.chains_to_design == None: #TODO this will likely break for sym mode

            if cfg.inference.symmetry!=None:#if we are doing symmetry:
                breaks = int(cfg.inference.symmetry[1:])
            else:
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

    # Check if we already have cfg.num_designs_rfdiff .pdb and *.trb files in cfg.rfdiff_out_path
    if len(glob.glob(os.path.join(cfg.rfdiff_out_path, "*.pdb"))) == cfg.num_designs_rfdiff:
        log.info(f"Found {cfg.num_designs_rfdiff} .pdb and .trb files in {cfg.rfdiff_out_path}. Skipping RFdiff.")
        if len(glob.glob(os.path.join(cfg.rfdiff_out_path, "*.trb"))) != cfg.num_designs_rfdiff:
            log.critical(f"Found {len(glob.glob(os.path.join(cfg.rfdiff_out_path, '*.trb')))} trb files in {cfg.rfdiff_out_path}, expected {cfg.num_designs_rfdiff}. \n Most likely, your previous simulation crashed while RfDiff was writing output files. Please, manually remove the .pdb file which does not have an associated .trb file in the {cfg.rfdiff_out_path} directory, then run the computation again.")
            raise Exception(f"Number of RfDiff .pdb files does not match RfDiff .trb files!") 
        #endif
        log.info("***************Skipping RFdiffusion altogether***************")
        return


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
            f'{cfg.pymol_python_path} {scripts_folder / "rechain.py"} "{pdb}" "{pdb}" --chain_break_cutoff_A {cfg.chain_break_cutoff_A}'
        )
        

    log.info("After rechaining")

def parse_additional_args(cfg, group):
    dodatniArgumenti = ""
    for k, v in (cfg.get(group, {}) or {}).items(): # or to allow for empty groups
        dodatniArgumenti += f" {k} {v}"
    return(dodatniArgumenti)

error_messages = ["blank", "Am Anfang", "Just before entering next cycle (before deleting cycle_dir)", "Deleted cycle_dir, then crashed. New wasn't created yet.", "New empty cycle_dir, but is_af_corrupted=1", "New empty cycle_dir and is_af_corrupted=0", "Empty cycle_dir, just after saving is_af_corrupted=0", "Now created a dir, but haven't saved is_af_corrupted to 0 yet.", "Created dir & set is_af_corrupted to 0.", "9: This one is nasty! Already moved from af to cycle_dir, but haven't saved iteration number yet. After restart, one file will be lost (overwritten). TODO: Use alternative checkpointing.", "10: Moved and saved iteration number.", "Removed MPNN dir", "After parse_multiple_cahins", "After proteinMPNN", "Saved is_af_corrupted to 1, but haven't removed the af dir yet", "15: AF is now corrupted and empty.", "One file written to /monomers. On restart, will it overwrite it with the same one?", "Uh-oh, it crashed just prior to running AF.", "18: Now it crashed just after running AF for one fasta file.", "19: Well, now it crashed after running AF for all fasta files. But before saving cycle checkpoint! Expected restart result: just re-run it, without copying/moving anything in the cycle_dir", "20: Crashed just before updating cycle checkpoint.", "21: Finished all cycles; before final_ops"]
def throw(id, cycle=0):
    log.info(f"Wanna crash? {id} ?= {crash_at_error} in cycle {cycle} ?= {crash_at_cycle}")
    if id == crash_at_error and cycle == crash_at_cycle:
        msg = error_messages[id]
        log.critical(f"Forced crash at {id}: {msg} in cycle {cycle}")
        raise Exception(f"Forced crash at {id}: {msg} in cycle {cycle}") 

def save_checkpoint(folder, piece, value):
    with open(os.path.join(folder, f"checkpoint_{piece}.txt"), "w") as f:
        f.write(str(value))
        f.flush()
    log.info(f"Saving checkpoint {piece} = {value}")

def get_checkpoint(folder, piece, default=0):
    # check if file exists
    if os.path.exists(os.path.join(folder, f"checkpoint_{piece}.txt")):
        with open(os.path.join(folder, f"checkpoint_{piece}.txt"), "r") as f:
            value = int(f.read() or 0)
            log.info(f"Reading checkpoint {piece}: {value}")
            return value # suppose only numbers (as strings) are there (user should not write to this file). In worst case, file is empty.
    else:
        log.info(f"Checkpoint {piece} doesn't exist. Returning default value {default}")
        return default

def do_cycling(cfg):
    """
    Runs a loop of (ProteinMPNN -> AF2 ->) `af2_mpnn_cycles` number of times.
    """
    log.info("Running do_cycling")
    start_cycle = get_checkpoint(cfg.output_dir, "cycle", 0)
    content_status = get_checkpoint(cfg.output_dir, "content_status", 0) # 0 ... fresh run; 1 ... can clear and copy; 2 ... can copy; 3 ... corrupted (partial new files in AF2)
    for cycle in range(start_cycle, cfg.af2_mpnn_cycles):
        print("cycleeeeee", cycle)
        dtimelog(f"Starting cycle {cycle}")

        trb_paths = None
        input_mpnn = cfg.rfdiff_out_dir # First cycle has this, other cycles overwrite this var
        throw(1, cycle)

        if not cycle == 0: # All other cycles get starting PDBs from AF2
            print(f"Cycle is not 0: {cycle}")
            cycle_directory = os.path.join(cfg.output_dir, "2_1_cycle_directory")

            if content_status == 1:
                throw(2, cycle)
                if os.path.exists(cycle_directory):
                    shutil.rmtree(cycle_directory) # We should not remove it on restart
                    throw(3, cycle)
                os.makedirs(cycle_directory, exist_ok=True)
                throw(4, cycle)
                save_checkpoint(cfg.output_dir, "content_status", 2) ## AF folder is ok to move its files to cycle_folder and cycle_folder is empty (ready for new files).
                content_status = 2
                throw(5, cycle)
            print(f"Nada: {cycle}")
            throw(20, cycle)

            if content_status == 2:
                af2_model_subdicts = glob.glob(os.path.join(cfg.af2_out_dir, "*"))
                log.info(f"AF2 model subdicts: {af2_model_subdicts}")

                #Access all af2 models and put them in one intermediate directory to get the same ored as in the 1st cycle (all pdbs in one directory)
                for model_subdict in af2_model_subdicts: 
                    log.info(f"Model subdict: {model_subdict}")
                    af2_pdbs = sorted(glob.glob(os.path.join(model_subdict, "T*.pdb")))
                    log.info(f"AF2 pdbs: {af2_pdbs}")
                    
                    rf_model_num = prosculpt.get_token_value(os.path.basename(model_subdict), "model_", "(\\d+)")
                    # Count number of files in cycle_directory which contain f"rf_{rf_model_num}"
                    number_of_copied_models = len([f for f in os.listdir(cycle_directory) if f"rf_{rf_model_num}" in f])
                    log.info(f"Number of already copied/present models: {number_of_copied_models}")
                    #last_iteration = max([int(prosculpt.get_token_value(f, "itr_", "(\\d+)")) for f in os.listdir(cycle_directory) if f"rf_{rf_model_num}" in f] or [0])

                    for i, af2_pdb in enumerate(af2_pdbs, start=number_of_copied_models):
                        af_model_num = prosculpt.get_token_value(os.path.basename(af2_pdb), "model_", "(\\d+)")
                        #if str(af_model_num) in args.af2_models:
                        # Rename pdbs to keep the model_num traceability with orginal rfdiff structure and enable filtering which models for next cycle

                        log.info(f"i: {i}, af2_pdb: {af2_pdb}, af_model_num: {af_model_num}, rf_model_num: {rf_model_num}")
                        shutil.move(af2_pdb, os.path.join(cycle_directory, f"rf_{rf_model_num}__model_{af_model_num}__cycle_{cycle}__itr_{i}__.pdb")) 
                                    #rf_ --> rfdiffusion structure number (in rfdiff outou dir)
                                    #model_ -> af2 model num, used for filtering which to cycle (preference for model 4)
                                    #itr_ -> to differentiate models in later cycles (5 pdbs for model 4 from rf 0 for example)
                                    # is it maybe possible to filter best ranked by af2 from the itr numbers?
                        throw(9, cycle)
                save_checkpoint(cfg.output_dir, "content_status", 3) ## AF folder is soon going to be in such a state that .pdb files should not be moved on restart (because you would have two different cycles in the same cycle_folder)
                content_status = 3

            input_mpnn = cycle_directory

            shutil.rmtree(cfg.mpnn_out_dir) # Remove MPNN dir so you can create new sequences
            throw(11, cycle)
            os.makedirs(cfg.mpnn_out_dir, exist_ok=True) # Create it again
            trb_paths = os.path.join(cfg.rfdiff_out_dir, "*.trb")
            print('trb_path is: ', trb_paths)
        #endif
        # All cycles run the same commands
        dtimelog(f"Cycle {cycle} helper scripts")
        run_and_log(
            f'{cfg.python_path_mpnn} {os.path.join(cfg.mpnn_installation_path, "helper_scripts", "parse_multiple_chains.py")} \
            --input_path={input_mpnn} \
            --output_path={cfg.path_for_parsed_chains}'       
        )
        throw(12, cycle)

        if cfg.chains_to_design:
            run_and_log(
                f"{cfg.python_path_mpnn} {os.path.join(cfg.mpnn_installation_path, 'helper_scripts', 'assign_fixed_chains.py')} \
                    --input_path={cfg.path_for_parsed_chains} \
                    --output_path={cfg.path_for_assigned_chains} \
                    --chain_list='{cfg.chains_to_design}'"
                    )
            
        fixed_pos_path = prosculpt.process_pdb_files(input_mpnn, cfg.mpnn_out_dir, cfg, trb_paths) # trb_paths is optional (default: None) and only used in non-first cycles
        # trb_paths is atm not used in process_pdb_files anyway -- a different approach is used (file.pdb -> withExtension .trb), which ensures the PDB and TRB files match.

        if cfg.inference.symmetry!=None:
            run_and_log(
                f'{cfg.python_path_mpnn} {os.path.join(cfg.mpnn_installation_path, "helper_scripts", "make_tied_positions_dict.py")} \
                --input_path={cfg.path_for_parsed_chains} \
                --output_path={cfg.path_for_tied_positions} \
                --homooligomer 1'    
            )         
            log.info(f"running symmetry")

        #_____________ RUN ProteinMPNN_____________
        # At first cycle, use num_seq_per_target from config. In subsequent cycles, set it to 1.
        dtimelog(f"Cycle {cycle} running ProteinMPNN")
        proteinMPNN_cmd_str = f'{cfg.python_path_mpnn} {os.path.join(cfg.mpnn_installation_path, "protein_mpnn_run.py")} \
            --jsonl_path {cfg.path_for_parsed_chains} \
            --fixed_positions_jsonl {cfg.path_for_fixed_positions} \
            {"--tied_positions_jsonl "+cfg.path_for_tied_positions if cfg.inference.symmetry!=None else ""} \
            --chain_id_jsonl {cfg.path_for_assigned_chains} \
            --out_folder {cfg.mpnn_out_dir} \
            --num_seq_per_target {cfg.num_seq_per_target_mpnn if cycle == 0 else 1} \
            --sampling_temp {cfg.sampling_temp} \
            --backbone_noise {cfg.backbone_noise} \
            --use_soluble_model  \
            --omit_AAs {cfg.omit_AAs} \
            {parse_additional_args(cfg, "pass_to_mpnn")} \
            --batch_size 1'
        
        run_and_log(proteinMPNN_cmd_str)
        throw(13, cycle)

        log.info("Preparing to empty af2 directory.")
        dtimelog(f"Cycle {cycle} preparation for af2")
        
        # af2 directory must be empty
        throw(14, cycle)
        shutil.rmtree(cfg.af2_out_dir)
        os.makedirs(cfg.af2_out_dir, exist_ok=True)
        throw(15, cycle)


        #________________ RUN AF2______________
        fasta_files = sorted(glob.glob(os.path.join(cfg.fasta_dir, "*.fa"))) # glob is not sorted by default
        print(fasta_files)

        monomers_fasta_dir = os.path.join(cfg.fasta_dir, 'monomers')
        if not os.path.exists(monomers_fasta_dir):
            os.makedirs(monomers_fasta_dir)
        
        #if symmetry - make fasta file with monomer sequence only
        for fasta_file in fasta_files:
            if cfg.inference.symmetry!=None or cfg.get("model_monomer", False):  #TODO: move this up one line because cfg symmetry is constant, not based on fasta_file
                sequences=[]
                for record in SeqIO.parse(fasta_file, "fasta"):
                    record.seq = record.seq[:record.seq.find('/')]  
                    idd = record.id            
                    record.id = 'monomer_'+idd      
                    descr = record.description
                    record.description = 'monomer_'+descr
                    sequences.append(record)
                SeqIO.write(sequences, os.path.join(os.path.dirname(fasta_file), "monomers", f"{os.path.basename(fasta_file)[:-3]}_monomer.fa"), "fasta") # It must be written to the correct subfolder already to prevent duplication on job restart (CHECKPOINTING): _monomer_monomer.fa
                print("File written to /monomers subfolder. "+fasta_file) #just for DEBUG 
                ## CHECKPOINTING: Right now, protMPNN is rerun after restart, which outputs different files into mpnn_out. The new files overwrite the old ones in /monomers, so it is always fresh.
                ## TODO: skip it (need more state variable checkpoints)
                throw(16, cycle)

        fasta_files = glob.glob(os.path.join(cfg.fasta_dir, "*.fa"))
        fasta_files += glob.glob(os.path.join(cfg.fasta_dir,'monomers', "*.fa"))
        fasta_files = sorted(fasta_files)
                
        print(fasta_files)


        for fasta_file in fasta_files: # Number of fasta files corresponds to number of rfdiff models
            print('RUNNING FASTA ',fasta_file)
            if cycle == 0:
                rf_model_num = prosculpt.get_token_value(os.path.basename(fasta_file), "_", "(\\d+)") #get 0 from _0.fa using reg exp
            else:
                rf_model_num = prosculpt.get_token_value(os.path.basename(fasta_file), "rf_", "(\\d+)") #get 0 from rf_0__model_1__cycle_2__itr_0__.pdb
                af2_model_num = prosculpt.get_token_value(os.path.basename(fasta_file), "model_", "(\\d+)") #get 1 from rf_0__model_1__cycle_2__itr_0__.pdb
            model_dir = os.path.join(cfg.af2_out_dir, f"model_{rf_model_num}") #create name for af2 directory name: model_0

            if 'monomer' in os.path.basename(fasta_file): #if we are doing symmetry - make an extra directory for modeling monomers with af2
                model_dir = os.path.join(model_dir,'monomers')
                
            prosculpt.change_sequence_in_fasta(cfg.rfdiff_pdb, fasta_file)
            if not os.path.exists(model_dir):
                os.makedirs(model_dir)

            model_order = str(cfg.model_order).split(",") # If only one number is passed, Hydra converts it to int
            num_models = len(model_order) #model order "1,2,3" -> num of models = 3
            
            
            if cfg.get("use_a3m", False): #If we're using a custom a3m, generate it for each sequence in the fasta
                print("Generating custom msa files")
                input_a3m_files=[]
                with open(fasta_file) as fasta_f:
                    a3m_filename=""
                    for line in fasta_f:       
                        if line[0]==">":
                            a3m_filename= ("".join([ c if (c.isalnum()  or c ==".") else "_" for c in line[1:-1] ]))+".a3m"
                        else:
                            mpnn_seq=line
                            trb_file=os.path.join(cfg.rfdiff_out_dir, "_"+str(rf_model_num)+".trb") #
                            print (a3m_filename)
                            custom_a3m_path=os.path.join(model_dir,a3m_filename)
                            prosculpt.make_alignment_file(trb_file,mpnn_seq,cfg.a3m_dir,custom_a3m_path)
                            input_a3m_files.append(custom_a3m_path)

                for a3m_file in input_a3m_files:
                    if cycle == 0: #have to run af2 differently in first cycle 
                        run_and_log(
                            f'source {cfg.af_setup_path} && {cfg.python_path_af2} {cfg.colabfold_setup_path} \
                                        --model-type alphafold2_multimer_v3 \
                                        --msa-mode mmseqs2_uniref_env \
                                        {a3m_file} {model_dir} \
                                        --model-order {cfg.model_order} \
                                        {parse_additional_args(cfg, "pass_to_af")} \
                                        --num-models {num_models}'
                                        ) #changed from single_sequence
                    else: 
                        for model_number in model_order:
                            if af2_model_num == model_number: # From the af2 model 4 you want only model 4 not also 2 and for 2 only 2 not 4 (--model_order "2,4")
                                num_models = 1
                                run_and_log(
                                    f'source {cfg.af_setup_path} && {cfg.python_path_af2} {cfg.colabfold_setup_path} \
                                        --model-type alphafold2_multimer_v3 \
                                        --msa-mode mmseqs2_uniref_env \
                                        {a3m_file} {model_dir} \
                                        --model-order {model_number} \
                                        {parse_additional_args(cfg, "pass_to_af")} \
                                        --num-models {num_models}'
                                        ) #Changed from single_sequence

            else:  #this is the normal mode of operations. Single sequence. 
                if cycle == 0: #have to run af2 differently in first cycle 
                    dtimelog(f"Cycle {cycle} running AF2")
                    throw(17, cycle)
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
                            dtimelog(f"Cycle {cycle}.{model_number} running AF2")
                            throw(17, cycle)
                            run_and_log(
                                f'source {cfg.af_setup_path} && {cfg.python_path_af2} {cfg.colabfold_setup_path} \
                                    --model-type alphafold2_multimer_v3 \
                                    --msa-mode single_sequence \
                                    {fasta_file} {model_dir} \
                                    --model-order {model_number} \
                                    {parse_additional_args(cfg, "pass_to_af")} \
                                    --num-models {num_models}'
                                    )#
            dtimelog(f"Cycle {cycle} finished AF2 for this fasta file {fasta_file}")
            throw(18, cycle)
            
        dtimelog(f"Cycle {cycle} finished AF2 for all fasta files")
        throw(19, cycle)
        save_checkpoint(cfg.output_dir, "content_status", 1) ## Content is ok to be copied to cycle_dir
        save_checkpoint(cfg.output_dir, "cycle", cycle+1) # Restart should start with next cycle.
        content_status = 1

        # msa single sequence makes sense for designed proteins
    throw(21)

def final_operations(cfg):
    log.info("Final operations")
    json_directories = glob.glob(os.path.join(cfg.af2_out_dir, "*"))

    print('do we already have csv files? If we need to re-run stats re have to delete them: ',os.listdir(cfg.output_dir) )
    if 'output.csv' in os.listdir(cfg.output_dir):
        os.remove(cfg.output_dir+'/output.csv')
    if 'scores_rg_charge_sap.csv' in os.listdir(cfg.output_dir):
        os.remove(cfg.output_dir+'/scores_rg_charge_sap.csv')
    if 'final_output.csv' in os.listdir(cfg.output_dir):
        os.remove(cfg.output_dir+'/final_output.csv')

    for model_i in json_directories:  # for model_i in [model_0, model_1, model_2 ,...]
        
        trb_num = prosculpt.get_token_value(os.path.basename(model_i), "model_", "(\\d+)") #get 0 from model_0 using reg exp

        if 'pdb_path' in cfg:
            prosculpt.rename_pdb_create_csv(cfg.output_dir, cfg.rfdiff_out_dir, trb_num, model_i, cfg.pdb_path, cfg.inference.symmetry, model_monomer=cfg.get("model_monomer", False))
        else:
            prosculpt.rename_pdb_create_csv(cfg.output_dir, cfg.rfdiff_out_dir, trb_num, model_i, control_structure_path=None, symmetry=cfg.inference.symmetry, model_monomer=cfg.get("model_monomer", False))
            
                
    csv_path = os.path.join(cfg.output_dir, "output.csv") #constructed path 'output.csv defined in rename_pdb_create_csv function
    run_and_log(
        f'{cfg.python_path} {scripts_folder / "scoring_rg_charge_sap.py"} {csv_path}'
        )
        
    scores_rg_path = os.path.join(cfg.output_dir, "scores_rg_charge_sap.csv") #'scores_rg_charge_sap.csv defined in scoring_rg_... script
    prosculpt.merge_csv(cfg.output_dir, csv_path, scores_rg_path) #, cfg.get("output_best", True), cfg.get("rmsd_threshold",3),cfg.get("plddt_threshold",90) used for filtering. not needed now.

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
    if 'pdb_path' in cfg:
        with open_dict(new_cfg["inference"]):
            new_cfg["inference"]["input_pdb"] = cfg.pdb_path

    log.info(f"Saving this config for RfDiff: {OmegaConf.to_yaml(new_cfg)}")
    # save to file
    with open(os.path.join(cfg.output_dir, "prosculpt2rfdiff.yaml"), "w") as f:
        f.write(OmegaConf.to_yaml(new_cfg))


    pass

import time
TIMEMEASURES = {}
TIMECALC = time.time()
PREVIOUS_MESSAGE = "Before app start"
def dtimelog(message, final=False):
    global TIMECALC
    global PREVIOUS_MESSAGE
    dt = round(time.time()-TIMECALC, 1)
    log.warning(f"* * * {PREVIOUS_MESSAGE} lasted {dt} s. Running {message} * * *")
    TIMEMEASURES[PREVIOUS_MESSAGE] = dt
    PREVIOUS_MESSAGE = message
    TIMECALC = time.time()
    
    if final:
        # Print all measures one per line
        for k,v in TIMEMEASURES.items():
            log.warning(f"{k} lasted {v} s.")

crash_at_error = 0 # Added through command line. Where should the test crash?
crash_at_cycle = 0 # Added through command line. In which cycle should the test crash?

@hydra.main(version_base=None, config_path="config", config_name="run")
def prosculptApp(cfg: DictConfig) -> None:
    log.info("Hydrić")
    log.info("The following configuration was passed: \n" + OmegaConf.to_yaml(cfg))

    #log.debug(f"Now in Hydra, cwd = {os.getcwd()}")

    dtimelog("general_config_prep")
    general_config_prep(cfg)
    config = HydraConfig.get()
    config_name = config.job.config_name
    log.info(f"config.runtime.config_sources: {config.runtime.config_sources}")
    config_path = [path["path"] for path in config.runtime.config_sources if path["schema"] == "file"][-1] # If run without config (through cmd), [1] is out of range. Assume passed config file is always the last one.
    shutil.copy(os.path.join(config_path,config_name+'.yaml'), os.path.join(cfg.output_dir, f"input.yaml"))

    global crash_at_error
    global crash_at_cycle
    crash_at_error = cfg.get("throw", -1)
    crash_at_cycle = cfg.get("crash_at_cycle", 0)

    if cfg.get('only_run_analysis',False):#only_run_analysis
        print('***Skip everything, go to final operations***')
        dtimelog("only final_operations")
        final_operations(cfg)

    else:
        if not cfg.get("skipRfDiff", False):
            dtimelog("pass_config_to_rfdiff")
            pass_config_to_rfdiff(cfg)
            dtimelog("run_rfdiff")
            run_rfdiff(cfg)
            dtimelog("rechain_rfdiff_pdbs")
            rechain_rfdiff_pdbs(cfg)
        else:
            log.info("*** Skipping RfDiff ***")
            dtimelog("skipping RfDiff")
            # Copy input PDB to RfDiff_output_dir and rename it to follow the token scheme
            shutil.copy(cfg.pdb_path, os.path.join(cfg.rfdiff_out_dir, "_0.pdb"))

        dtimelog("do_cycling")
        do_cycling(cfg)
        dtimelog("final_operations")
        final_operations(cfg)
        dtimelog("Finished", True)


if __name__ == "__main__":
    print("File: ", __file__)
    a = pathlib.Path(__file__).resolve().parent
    log.info(f"Before running MYapp, cwd = {os.getcwd()}")
    prosculptApp()