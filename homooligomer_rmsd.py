import os
from Bio.PDB import PDBParser, PDBIO, Superimposer, PPBuilder
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.PDB import CEAligner
from Bio.PDB import MMCIFParser
import scipy
import numpy as np
from scipy.spatial import distance
import pandas as pd
import argparse

#we want to align corresponding chains in mobile and target structures
#to do that, we see which chain of mobile structure gives the lowest rmsd in alignment with first chain of target
#then we align these chains
#after that, we calculate distances between centers of mass of all the chains to find corresponding pairs (the closest)
#then we rearrange the chains of mobile structure in right order and align two structures

def align_chain_A(path_to_target, path_to_mobile, parser): 
    structure_target = parser.get_structure("target", path_to_target)
    structure_mobile = parser.get_structure("mobile", path_to_mobile)
    mobile_chains = list(structure_mobile.get_chains())
    target_chains = list(structure_target.get_chains())

    target_chain_res = list(target_chains[0].get_residues())
    target_chain_res = [ind['CA'] for ind in target_chain_res] 

    #let's find which chain in mobile gives the smallest rmsd in alignment with chain[0] of target
    list_rmsd_chains = []
    for n in mobile_chains:
        mobile_chain_res = list(n.get_residues())
        mobile_chain_res = [ind['CA'] for ind in mobile_chain_res]

        superimposer = Superimposer()
        superimposer.set_atoms(target_chain_res, mobile_chain_res)
        superimposer.apply(structure_mobile.get_atoms())
        list_rmsd_chains.append(superimposer.rms)

    index_min = np.argmin(list_rmsd_chains)

    #align chains with the best correspondence
    mobile_chain_res = list(mobile_chains[index_min].get_residues())
    mobile_chain_res = [ind['CA'] for ind in mobile_chain_res]

    superimposer = Superimposer()
    superimposer.set_atoms(target_chain_res, mobile_chain_res)
    superimposer.apply(structure_mobile.get_atoms())

    return structure_target, structure_mobile


#calculate center of mass of a chain
def chains_com_coords(chainlist):
    chains_dict = {}
    for i in range(len(chainlist)):
        chainname = chainlist[i].get_id()
        chain_com = chainlist[i].center_of_mass()
        chains_dict[chainname]=chain_com
    return chains_dict


def align_oligomers(path_to_target, path_to_mobile):
    
    parser = PDBParser(PERMISSIVE=1)
    structure_target, structure_mobile = align_chain_A(path_to_target, path_to_mobile, parser)
    target_chains = list(structure_target.get_chains())
    mobile_chains = list(structure_mobile.get_chains())

    target_chains_com = chains_com_coords(target_chains)
    mobile_chains_com = chains_com_coords(mobile_chains)

    distances_array = []
    for key_rfdiff in target_chains_com:
        distances_array.append([])
        for key_af in mobile_chains_com:
            dst = distance.euclidean(target_chains_com[key_rfdiff], mobile_chains_com[key_af])
            distances_array[-1].append(dst)

    #get dictionaries of corresponding chain names 
    corresponding = np.argmin(distances_array, axis=1)
    rename_dict_1 = {}
    rename_dict_2 = {}
    for i, j in zip(corresponding, range(len(corresponding))):
        str1=list(target_chains_com.keys())
        str2=list(mobile_chains_com.keys())
        str3=[n.lower() for n in str2]
        rename_dict_1[str1[i]]=str3[j]
        rename_dict_2[str3[j]]=str2[j]

    #rename chains of mobile structure
    #if the structures are very different, multiple chains in mobile can be close to target 
    #in that case, return rmsd = -1
    try:
        for model in structure_mobile:
            for chain in model:
                old_name = chain.get_id()
                chain.id = rename_dict_1[old_name]
        for model in structure_mobile:
            for chain in model:
                old_name = chain.get_id()
                chain.id = rename_dict_2[old_name] 
    except KeyError:
        rmsd = -1
        return rmsd
    
    #get a list of residues for alignment for mobile structure
    target_res_total = list(structure_target.get_residues())
    mobile_res_total = []
    for i in list(structure_target.get_chains()):
        for k in mobile_chains:
            if k.id == i.id:
                mobile_res_total += list(k.get_residues())   
    target_res_total = [ind['CA'] for ind in target_res_total]
    mobile_res_total = [ind['CA'] for ind in mobile_res_total]

    superimposer = Superimposer()
    superimposer.set_atoms(target_res_total, mobile_res_total)
    superimposer.apply(structure_mobile.get_atoms())
    rmsd = superimposer.rms

    return rmsd


if __name__ == '__main__':
    # Define the command-line arguments
    argparser = argparse.ArgumentParser(description='calculate rmsd of two oligomers')
    argparser.add_argument('--path_to_target', '-t', type=str,
                        help='path to the pdb of target protein (the one to align to)')
    argparser.add_argument('--path_to_mobile', '-m',  type=str,
                        help='path to the pdb of mobile protein')
    # Parse the arguments
    args = argparser.parse_args()
    path_to_target = args.path_to_target
    path_to_mobile = args.path_to_mobile
    
    # Call the align_oligomers function with the parsed arguments
    rmsd = align_oligomers(path_to_target, path_to_mobile)
    print('rmsd: ',rmsd)
