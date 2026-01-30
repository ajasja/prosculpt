from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from Bio.PDB import *
import numpy as np
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import argparse
import os

# Initialize PyRosetta
init()

# Define functions for Rg, charge, and SAP score calculation
def calculate_rg(chain):
    """
    Calculates the radius of gyration (Rg) of a protein chain using only the alpha carbons (CA).
    
    Parameters:
        chain (Bio.PDB.Chain.Chain): The protein chain to calculate the Rg of.
        
    Returns:
        The Rg of the protein chain.
    """
    # Get the coordinates of all alpha carbons in the chain
    ca_atoms = [atom for atom in chain.get_atoms() if atom.get_name() == 'CA']
    ca_coords = np.array([atom.get_coord() for atom in ca_atoms])

    # Calculate the Rg of the protein
    Rg = np.sqrt(np.mean(np.sum((ca_coords - np.mean(ca_coords, axis=0))**2, axis=-1)) + 1e-8)

    return Rg

def calculate_charge(chain, ph=7.4):
    """
    Calculates the charge of the protein chain at a specific pH.
    
    Parameters:
        chain (Bio.PDB.Chain.Chain): The protein chain to calculate the charge of.
        ph (float): The pH value to calculate the charge at. Default is 7.4.
        
    Returns:
        The charge of the protein at the specified pH.
    """
    # Extract the sequence of amino acid residues in the chain
    sequence = ''
    for residue in chain.get_residues():
        resname = residue.get_resname()
        if resname in aa3:
            sequence += aa1[aa3.index(resname)]
        else:
            print(f"Skipping residue {resname} because it is not a standard amino acid.")
            continue

    # Create a ProteinAnalysis object from the sequence
    protein_analysis = ProteinAnalysis(sequence)

    # Calculate the charge of the protein at a specific pH
    charge = protein_analysis.charge_at_pH(ph)

    return charge

def calculate_sap_score(pdb_file_path, chain="B"):
    
    # Load the PDB file into a pose object
    pose = pose_from_pdb(pdb_file_path)

    # Select only chain B using a SwitchChainOrder mover
    select_chain = XmlObjects.static_get_mover(f'<SwitchChainOrder name="so" chain_order="{chain}"/>')
    chain = pose.clone()
    select_chain.apply(chain)

    # Calculate the SAP score for chain B
    sap_score_metric = XmlObjects.static_get_simple_metric('<SapScoreMetric name="sap_metric"/>')
    sap_score_value = sap_score_metric.calculate(chain)

    # Return the SAP score value
    #sap_score_value = sap_score_metric.get(1)
    return sap_score_value

# Set up command line arguments
parser = argparse.ArgumentParser(description='Calculate Rg, charge, and SAP score for a list of PDB models.')
parser.add_argument('input_csv', type=str, help='Path to input CSV file')
parser.add_argument('--output_dir', type=str, help='Path to output directory. If not provided, output files will be written to the same directory as the input CSV.')
#add a flag minimize structure and overwrite the pdb. minimization before metrics
# Parse command line arguments
args = parser.parse_args()

# Read input CSV file
data = pd.read_csv(args.input_csv)

# Set output directory
if args.output_dir:
    output = args.output_dir
else:
    output = os.path.dirname(args.input_csv)

"""
# Calculate Rg, charge, and SAP score for each PDB model and write results to output CSV file
with open(f"{output}/rosetta_scores.csv", "a") as f:
    f.write("model_path,rg,charge,sap\n")
    for pdb in data["model_path"]:
        parser = PDBParser()
        structure = parser.get_structure('protein', pdb)
        chain = structure[0]['A']
        rg = calculate_rg(chain)
        charge = calculate_charge(chain, ph=7.4)
        sap = calculate_sap_score(pdb, "A")
        f.write(f"{pdb},{rg},{charge},{sap}\n")
"""
## Calculate Rg, charge, and SAP score for each PDB model and write results to output CSV file
for pdb in data['model_path']:
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb)
    chain = structure[0]['A']
    dictionary = {'model_path' : pdb,
                'rg': calculate_rg(chain),
                'charge': calculate_charge(chain, ph=7.4),
                'sap': calculate_sap_score(pdb, "A")}
    df = pd.DataFrame(dictionary, index=[0])

    path_csv = os.path.join(output, "rosetta_scores.csv")
    df.to_csv(path_csv, mode='a', header=not os.path.exists(path_csv), index=False) # header=not os.path.exists(path_csv) determines if including column names (header) in the CSV
