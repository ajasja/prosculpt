
import pickle
from Bio.PDB import PDBParser, PDBIO, Superimposer, PPBuilder
from Bio import SeqIO
from Bio.Seq import Seq
import json
import glob
import os
import numpy as np
import pandas as pd
import shutil
from pathlib import Path


def calculate_RMSD_linker_len (trb_path, af2_pdb, starting_pdb):        
    # First calculate RMSD between input protein and AF2 generated protein
    # Second calcualte number of total generated AA by RFDIFF 
    #   - if designing only in one location the number is equal linker length      
    
        with open(trb_path, 'rb') as f:
                trb_dict = pickle.load(f)
        parser = PDBParser(PERMISSIVE = 1)
        structure_control = parser.get_structure("control", starting_pdb)
        structure_designed = parser.get_structure("designed", af2_pdb)

        # Get different data from trb file depending on the fact if designing a monomer (one chain) or heteromer
        # complex_con_rex_idx present only if designing heteromer
        # Data structure: con_ref_idx0 = [0, 1, 3, 3, ...] 
        #   Info: array of input pdb AA indices starting 0 (con_ref_pdb_idx), and where they are in the output pdb (con_hal_pdb_idx)
        #   In complex_con_hal_idx0 there is no chain info however RFDIFF changes pdb indeces to go from 1 to n (e.g. 1st AA in chain B has idx 34)
        if 'complex_con_ref_idx0' in trb_dict:
            residue_data_control = trb_dict['complex_con_ref_idx0']
            residue_data_designed = trb_dict['complex_con_hal_idx0']
        else:
            residue_data_control = trb_dict['con_ref_idx0']
            residue_data_designed = trb_dict['con_hal_idx0']

        control_res = list(structure_control.get_residues()) #obtain a list of all the residues in the structure, structure_control is object
        control_res = [control_res[ind]['CA'] for ind in residue_data_control] #retrieve the residue with the corresponding index from control_res
        #printed_ref = [res.get_resname() for res in control_res]

        designed_res = list(structure_designed.get_residues())
        designed_res = [designed_res[ind]['CA'] for ind in residue_data_designed]
        
        #printed_desi = [res.get_resname() for res in designed_res]

        #print(list(x == y for x, y in zip(printed_ref, printed_desi)))
        
        trb_help = list(trb_dict['inpaint_str'])
        linker_indeces = [boolean for boolean in trb_help if boolean == False] #calculate linker length here - convenient
        linker_length = len(linker_indeces)


        superimposer = Superimposer()
        superimposer.set_atoms(control_res, designed_res)
        

        superimposer.apply(structure_designed.get_atoms())
        rmsd = superimposer.rms
        io=PDBIO()
        io.set_structure(structure_designed)
        io.save("af2_pdb_2.pdb")

        return (round(rmsd, 1), linker_length)



def get_token_value(astr, token, regular_expression): #"(\d*\.\d+|\d+\.?\d*)"
    """returns value next to token"""
    import re
    regexp = re.compile(f'{token}{regular_expression}')
    match = regexp.search(astr)
    #if match == None:
     #   match = "/"
    return match.group(1)


def merge_csv(output_dir, output_csv, scores_csv):
    # read csv files
    df1 = pd.read_csv(scores_csv)
    df2 = pd.read_csv(output_csv)

    # merge dataframes on 'model_path' column
    merged_df = pd.merge(df1, df2, on='model_path')
    # drop duplicate 'model_path' column (if it exists)
    merged_df = merged_df.loc[:, ~merged_df.columns.duplicated()]

    # save merged dataframe to csv file
    merged_df.to_csv(f'{os.path.join(output_dir, "final_output.csv")}', index=False)

def rename_pdb_create_csv(output_dir, rfdiff_out_dir, trb_num, model_i, control_structure_path):

    # Preparing paths to acces correct files
    model_i = os.path.join(model_i, "") # add / to path to access json files within
    dir_renamed_pdb = os.path.join(os.path.dirname(output_dir), "final_pdbs")
    os.makedirs(dir_renamed_pdb, exist_ok=True) # directory is created even if some or all of the intermediate directories in the path do not exist

    trb_file = os.path.join(rfdiff_out_dir, f"_{trb_num}.trb") #name of corresponding trb file 
    rfdiff_pdb_path = os.path.join(rfdiff_out_dir, f"_{trb_num}.pdb")

    json_files = glob.glob(os.path.join(model_i, 'T*000.json'))

    for json_file in json_files: #in af2 model_i directory for T_...json file in [all T_...json files] 
        # This is done for each model_i directory therefore for each rfdiff pdb 
        # There are 5 T_...jsons per 1 mpnn seq of the rfdiff model

        with open(json_file, 'r') as f:
            params = json.load(f)

        # Handle filenames correctly to get the T_...pdb that corresponds to the T_...json
        json_filename = os.path.basename(json_file)
        json_dirname = os.path.dirname(json_file)
        json_newname = os.path.join(json_dirname, json_filename.replace("scores", "unrelaxed"))
        model_pdb_file = os.path.splitext(json_newname)[0] + '.pdb' #T_0.1__sample_2__score_0.5830__global_score_0.8339__seq_recovery_02000_unrelaxed_rank_004_alphafold2_multimer_v3_model_3_seed_000.pdb
        
        

        #Extract relevant data. Files used: json file of specific af2 model, specific af2 pdb,  trb file of rfdiff model (1 for all AF2 models from same rfdiff pdb)  
        plddt_list = params['plddt']
        plddt = int(np.mean(plddt_list))
        rmsd, linker_length	= calculate_RMSD_linker_len (trb_file, model_pdb_file, control_structure_path)
        pae = round((np.mean(params['pae'])), 2)
        
        #tracebility
        output_num = os.path.basename(output_dir)
        af2_model =  get_token_value(json_filename, '_model_', "(\d*\.\d+|\d+\.?\d*)")

        # Create a new name an copy te af2 model under that name into the output directory
        new_pdb_file = f"link_{linker_length}__plddt_{plddt}__rmsd_{rmsd}__pae_{pae}__out_{output_num}__rf_{trb_num}__af_model_{af2_model}_.pdb"
            #out -> 00 -> number of task
            #rf -> 01 -> number of corresponding rf difff model
            #af_model -> 4 -> number of the af model (1-5), can be set using --model_order flag 
        new_pdb_path = os.path.join(dir_renamed_pdb, new_pdb_file)

        try:
            shutil.copy2(model_pdb_file, new_pdb_path)
        except OSError as e:
            print(f"Error copying {model_pdb_file} to {new_pdb_file}: {e}")

        p = PDBParser()
        structure = p.get_structure("model_seq", new_pdb_path)
        ppb = PPBuilder()
        
        seq = ""
        for pp in ppb.build_peptides(structure):
            seq += f":{pp.get_sequence().__str__()}"



        dictionary = {'link_lenght': get_token_value(new_pdb_file, 'link_', "(\d*\.\d+|\d+\.?\d*)" ),
                'plddt': get_token_value(new_pdb_file, '__plddt_', "(\d*\.\d+|\d+\.?\d*)"),
                'RMSD': get_token_value(new_pdb_file, '__rmsd_', "(\d*\.\d+|\d+\.?\d*)"),
                'pae': get_token_value(new_pdb_file, '__pae_', "(\d*\.\d+|\d+\.?\d*)"),
                'model_path': new_pdb_path,
                'sequence' : seq[1:],
                'af2_json' : json_file,
                'af2_pdb' : model_pdb_file,
                'path_rfdiff': rfdiff_pdb_path }  #MODEL PATH for scoring_rg_... #jsonfilename for traceability

        df = pd.DataFrame(dictionary, index=[0])
        path_csv = os.path.join(output_dir, "output.csv")
        df.to_csv(path_csv, mode='a', header=not os.path.exists(path_csv), index=False)


# vsaka tabelca za svoj task
# funkcija na koncu, ki vse združi
# sestavljanje pathov je ok, dodaj nov column s pathom do rfdif
        
        
        




def create_dataframe(path_to_files, output_dir): # path = r'content/*partial.pdb' 
    #takes path  to renamed pdbs 
    all_files = glob.glob(os.path.join(path_to_files, "*.pdb"))
    list_of_dicts = []

    for file_name in all_files:
        p = PDBParser()
        structure = p.get_structure("model_seq", file_name)
        ppb = PPBuilder()
        
        seq = ""
        for pp in ppb.build_peptides(structure):
            seq += f":{pp.get_sequence().__str__()}"
        print(seq)
        dictionary = {'link_lenght': get_token_value(file_name, 'link_', "(\d*\.\d+|\d+\.?\d*)" ),
                'plddt': get_token_value(file_name, '__plddt_', "(\d*\.\d+|\d+\.?\d*)"),
                'RMSD': get_token_value(file_name, '__rmsd_', "(\d*\.\d+|\d+\.?\d*)"),
                'pae': get_token_value(file_name, '__pae_', "(\d*\.\d+|\d+\.?\d*)"),
                'model_path': file_name,
                'sequence' : seq,
                'rfdiff model': get_token_value(file_name, '__rf_', "(\d*\.\d+|\d+\.?\d*)")}  #MODEL PATH for scoring_rg_... #jsonfilename for traceability
        
        list_of_dicts.append(dictionary)

    #columns = ['link_length', 'plddt', 'loop_plddt', 'RMSD', 'model_path', 'sequence', 'score_traceb']
    df = pd.DataFrame(list_of_dicts)
    path_csv = os.path.join(os.path.dirname(output_dir), "output.csv")
    df.to_csv(path_csv, mode='a', header=not os.path.exists(path_csv), index=False)

class NumpyInt64Encoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.int64):
            return int(obj)
        return super(NumpyInt64Encoder, self).default(obj)


def process_pdb_files(pdb_path: str, out_path: str, trb_paths = None):
    fixpos = {}
    pdb_files = Path(pdb_path).glob("*.pdb")

    for pdb_file in pdb_files:
        
        pdb_basename = pdb_file.stem
        trb_file = pdb_file.with_suffix(".trb")

        if not trb_file.exists():
            rf_model_num = get_token_value(os.path.basename(pdb_file), "rf_", "(\d+)")
            trb_file = os.path.join(os.path.dirname(pdb_file), f"_{rf_model_num}.trb")
            trb_file = trb_file.replace("2_1_cycle_directory", "1_rfdiff") 
            print(f"TRB file not found for {pdb_basename}. CAUTION, using composed path")
        
        print(pdb_file, trb_file)

        # We need to renumber fixed resids: each chain should start with 1 
        chainResidOffset = {}
        parser = PDBParser()
        chainsInPdb = parser.get_structure("protein", pdb_file).get_chains()
        for let_chain in chainsInPdb:
            chainLetter = let_chain.get_id()
            startingNo = int(next(let_chain.get_residues()).get_id()[1])
            chainResidOffset.setdefault(chainLetter, startingNo-1)
                    

        with open(trb_file, 'rb') as f:
            trb_data = pickle.load(f)

        # fixed_res should never be empty, otherwise ProteinMPNN will throw a KeyError fixed_position_dict[b['name']][letter]
        # We need to set blank fixed_res for each generated chain (based on contig).
        contig = trb_data["config"]["contigmap"]["contigs"][0]
        abeceda = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        breaks = contig.count("/0 ") + 1
        fixed_res = dict(zip(abeceda, [[]]*breaks))
        print(f"Fixed res (according to contig chain breaks): {fixed_res}")
        
        if 'complex_con_hal_pdb_idx' in trb_data:
            con_hal_idx = trb_data.get('complex_con_hal_pdb_idx', []) #con_hal_pdb_idx #complex_con_hal_pdb_idx
        else:
            con_hal_idx = trb_data.get('con_hal_pdb_idx', [])
        # Process con_hal_idx to extract chain ids and indices
        
        # This is only good if multiple chains due to symmetry: all of them are equal; ProteinMPNN expects fixed_res as 1-based, resetting for each chain.
        for chain, idx in con_hal_idx:
            # If there are multiple chains, reset the auto_incrementing numbers to 1 for each chain (subtract offset)
            fixed_res.setdefault(chain, []).append(idx - chainResidOffset[chain]) 
            # RfDiff outputs multiple chains if contig has /0 (chain break)
            
        fixpos[pdb_basename] = fixed_res
    
    #print("_________trb data____", trb_data)
    
    
    # print("_________ fix pos_________", fixpos)
    file_path = os.path.join(out_path, "fixed_pdbs.jsonl")
    # Save the fixpos dict as a JSON file
    with open(file_path, "w") as outfile:
        json.dump(fixpos, outfile, cls=NumpyInt64Encoder)


    return file_path

def get_chains_seq(pdb_file):

    """
    Extract sequences of all chains in a protein from a PDB file.
    Args:
        pdb_file (str): The path to the PDB file.

    Returns:
        list: A list of strings of seqenceces of all chains.

        Assumptions:
        - RfDIFF forms the new pdb in a way that the connected helices are now chain A, regardless of how they were identified before
        - Other chains are now B, C,... regardles if before they were A

    """
    # Parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)

    # Extract polypeptides and their sequences
    pp_builder = PPBuilder()
    other_chains_sequences = []

    for pp in pp_builder.build_peptides(structure):
        sequence = pp.get_sequence()
        seqs= Seq(sequence)
        other_chains_sequences.append(seqs)

    return other_chains_sequences

def read_fasta_file(fasta_file):
    """
    Read a FASTA file using Biopython.

    Args:
        fasta_file (str): The path to the FASTA file.

    Returns:
        list: A list of SeqRecord objects containing the sequences from the FASTA file.
    """
    sequences = []

    with open(fasta_file, "r") as file_handle:
        for record in SeqIO.parse(file_handle, "fasta"):
            sequences.append(record)
        #sequences = list(SeqIO.parse(file_handle, "fasta")) # better, for loop ni potreben
    return sequences

def change_sequence_in_fasta (pdb_file, mpnn_fasta):
        """
            #function adds : for chainbraks and chain sequences of all other chains
        sequences_all_chains = get_chains_seq(pdb_file)
        sequences_other_chains = sequences_all_chains[1:] 
        #sequences of other chains because mpnn fasta has only chain A and no other chain seqs

        sequences_mpnn = read_fasta_file(mpnn_fasta) 
        #funkcija za drop duplicates, napiši funkcijo
        for seq_record in sequences_mpnn:
                new_sequence = seq_record.seq
                for other_chain in sequences_other_chains:
                        new_sequence += f":{other_chain}"
                seq_record.seq = Seq(new_sequence)

        with open(mpnn_fasta, "w") as output:
                SeqIO.write(sequences_mpnn, output, "fasta")
        """
        # Read text file, then replace all slashes with colons
        with open(mpnn_fasta, "r+") as f:
            text = f.read()
            text = text.replace("/", ":")
            f.seek(0)
            f.write(text)
            f.truncate()
        # I doubt we need to manually add "other" chains -- according to tutorials/examples, contig should always specify all chains (even if unchanged): [A1-30/6-6/B1-30/0 C1-100/0 D1-25/0 E1-123]


def match_linker_length(trb_path):
    import re

    with open(trb_path, 'rb') as f:
        trb_dict = pickle.load(f)

    input_string = ' '.join(trb_dict['sampled_mask']) #' '.join(['A1-30/6-6/C1-30/0', 'D1-30/0', 'B1-30/0'])

    pattern = r'(?<=/)(\d+)(?=-\d+/)'
    """
        (?<=/): Positive lookbehind assertion -> sub-pattern is preceded by a '/'.
        (\d+): Capture group that matches one or more digits.
        (?=-\d+/): Positive lookahead assertion -> sub-pattern is followed by a '-' and one or more digits, then a '/'.
    """
    match = re.search(pattern, input_string)
    if match:
        link = match.group(1)
        return link  # Output: 6
    else:
        print("No match found in linker length function")

def find_correct_trb(model_i):
    json_files = glob.glob(os.path.join(model_i, '_*000.json'))
    
    print(json_files)
    for json_file in json_files:
        json_filename = os.path.basename(json_file)
        number_trb = get_token_value(json_filename, "_", "(\d+)(?=_)") # stop at first digit -> (?=_)
        trb_file = os.path.join(f"_{number_trb}.trb")
    return trb_file