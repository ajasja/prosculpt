import os
from Bio.PDB import PDBParser, PDBIO, Superimposer, Structure
from scipy.spatial import distance
import numpy as np
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
        # chainname = chainlist[i].get_id()
        chain_com = chainlist[i].center_of_mass()
        chains_dict[i]=chain_com
    return chains_dict

#if set to true, saves aligned structures to a file
def save_structures(structure_mobile, path_to_target, path_to_mobile, rmsd, output_pdb):#dir->pdb; 
    if output_pdb:
        pdb_file = output_pdb
    else:
        pdb_file = f"{os.path.basename(path_to_target)}_{os.path.basename(path_to_mobile)}_rmsd_{rmsd}.pdb"
    io=PDBIO()
    io.set_structure(structure_mobile)
    io.save(pdb_file)

#if set to true, aligns a monomer to every chain of oligomer and returns the best rmsd
def align_monomer(path_to_target, path_to_mobile, save_aligned,  output_pdb=None):

    parser = PDBParser(PERMISSIVE=1)
    structure_target = parser.get_structure("target", path_to_target)
    structure_mobile = parser.get_structure("mobile", path_to_mobile)
    
    target_chains = list(structure_target.get_chains())
    mobile_chain_res = list(structure_mobile.get_residues())
    mobile_chain_res = [ind['CA'] for ind in mobile_chain_res]

    list_rmsd_chains = []
    for n in target_chains:
        target_chain_res = list(n.get_residues())
        target_chain_res = [ind['CA'] for ind in target_chain_res]

        superimposer = Superimposer()
        superimposer.set_atoms(target_chain_res, mobile_chain_res)
        superimposer.apply(structure_mobile.get_atoms())
        list_rmsd_chains.append(superimposer.rms)
    rmsd = np.min(list_rmsd_chains)
    print('rmsd: ', rmsd)
    
    #align chains with the best correspondence (only aiign again so that we save the best-aligned option)
    if save_aligned==True:
        index_min = np.argmin(list_rmsd_chains)
        target_chain_res = list(target_chains[index_min].get_residues())
        target_chain_res = [ind['CA'] for ind in target_chain_res]
        superimposer = Superimposer()
        superimposer.set_atoms(target_chain_res, mobile_chain_res)
        superimposer.apply(structure_mobile.get_atoms())
        rmsd = superimposer.rms
        save_structures(structure_mobile, path_to_target, path_to_mobile, rmsd,output_pdb)
    
    return rmsd


def align_oligomers(path_to_target, path_to_mobile, save_aligned, output_pdb=None):
    print('output_pdb ',output_pdb)
    print('path_to_target ',path_to_target)
    
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

    #get dictionary of corresponding chain indexes 
    corresponding = np.argmin(distances_array, axis=1)

    target_res_total = list(structure_target.get_residues())
    mobile_res_total = []
    for i in range(len(list(structure_target.get_chains()))):
        k = list(structure_mobile.get_chains())[corresponding[i]]
        mobile_res_total += list(k.get_residues())
    target_res_total = [ind['CA'] for ind in target_res_total]
    mobile_res_total = [ind['CA'] for ind in mobile_res_total]

    superimposer = Superimposer()
    superimposer.set_atoms(target_res_total, mobile_res_total)
    superimposer.apply(structure_mobile.get_atoms())
    rmsd = superimposer.rms
    print('rmsd: ', rmsd)
    if save_aligned==True:
        save_structures(structure_mobile, path_to_target, path_to_mobile, rmsd, output_pdb)

    return rmsd


if __name__ == '__main__':
    # Define the command-line arguments
    argparser = argparse.ArgumentParser(description='calculate rmsd of two oligomers')
    argparser.add_argument('--path_to_target', '-t', type=str,
                        help='path to the pdb of target protein (the one to align to)')
    argparser.add_argument('--path_to_mobile', '-m',  type=str,
                        help='path to the pdb of mobile protein')
    argparser.add_argument('--output_pdb', '-out', type=str, default=None,
                        help='directory/filename - where to save a file with aligned mobile structure') #if we pass, then saves there; if we do not, then saves a combination of names of target and mobile
    argparser.add_argument('--save_aligned', '-s',  type=bool, default=True, 
                        help='if you want to save the aligned mobile structure')
    argparser.add_argument('--align_monomer', '-am',  type=bool, default=False,
                        help='do you want to align a monomer?')
    
    # Parse the arguments
    args = argparser.parse_args()
    path_to_target = args.path_to_target
    path_to_mobile = args.path_to_mobile
    output_pdb = args.output_pdb
    save_aligned = args.save_aligned
    ismonomer = args.align_monomer
    
    # Call the align_oligomers function with the parsed arguments
    if ismonomer:
        rmsd = align_monomer(path_to_target, path_to_mobile,output_pdb, save_aligned)
    else:
        rmsd = align_oligomers(path_to_target, path_to_mobile,output_pdb, save_aligned)
    

