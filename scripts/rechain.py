"""Pymol script to assign chain IDs based on chain breaks. Made by ajasja.ljubetic@gmail.com"""

try:
    import pymol
    from pymol import cmd
    from pymol import stored
    from pymol import CmdException
except:
    cmd = None
    pymol = None
    print("Pymol is most likely not installed")
import os
import argparse

#/home/aljubetic/conda/envs/pyro/bin/python


def get_pymol_name(file_name):
    """Gets pymol name from file name"""
    name = os.path.basename(file_name)
    name = os.path.splitext(name)[0]
    #make it work with pdb.gz
    if name.endswith('.pdb'):
        name = name[:-4]
    name = cmd.get_legal_name(name)
    return name

def get_chainbreaks(objname_or_file, cutoff_A=2, cmd=None, delete_all_before_load=False):
    """Returns the resi where the chain break is. The next residue is in a different place 
    N --chain_break-- N+1  .Returns N."""
    import numpy as np
    if cmd is None:
        import pymol
        cmd = pymol.cmd

    if os.path.exists(objname_or_file):
        objname = get_pymol_name(objname_or_file)
        if delete_all_before_load:
            cmd.delete('all')
        cmd.load(objname_or_file, object=objname)
    else:
        objname = objname_or_file 

    C_atoms = cmd.get_coords(objname+" and name C", 1)
    N_atoms = cmd.get_coords(objname+" and name N", 1)
    print(objname)
    print(C_atoms[:-1])
    print(N_atoms[1:])
    #subastract the C from the N of the next residue
    distances = np.sum((C_atoms[:-1]-N_atoms[1:])**2, axis=1)
    #len(distances), min(distances), np.mean(distances) ,max(distances)
    breaks = distances > cutoff_A**2
    return breaks.nonzero()[0]+1

def re_split_chains(in_name, out_name, chain_break_cutoff_A=2, cmd=None, center_outer_segid=None):
    """Splits chains by chain-break. Residues must be numbered 1-N?"""
    if cmd is None:
        import pymol
        cmd = pymol.cmd

    chain_breaks = get_chainbreaks(in_name, cutoff_A=chain_break_cutoff_A, cmd=cmd, delete_all_before_load=True)
    chain_breaks = [0] + list(chain_breaks) +[' '] #First chain is at position 0
    print(chain_breaks)
    chainsID = "ABCDEFGHIJKLMNOPRSTUVXYZ0123456789abcdefghijklmnoprstuvxyz!@#$%^&*()_+"
    for n in range(len(chain_breaks)-1):
        #Tule je malo psevdokoda manjka za n=0 in n=len(cb)
        cmd.alter(f"resi {chain_breaks[n]+1}-{chain_breaks[n+1]}", f"chain='{chainsID[n]}'") #####_______Explain better____
        print(f"fresi {chain_breaks[n]+1}-{chain_breaks[n+1]}")
     
    cmd.do(f"save {out_name}")
    return out_name

if __name__ == '__main__':
    # Define the command-line arguments
    parser = argparse.ArgumentParser(description='Split chains by chain-break.')
    parser.add_argument('in_name', help='Input file name')
    parser.add_argument('out_name', help='Output file name')
    parser.add_argument('--chain_break_cutoff_A', type=float, default=2, help='Chain break cutoff in angstroms (default: 2)')

    
    # Parse the arguments
    args = parser.parse_args()
    
    # Call the re_split_chains function with the parsed arguments
    re_split_chains(args.in_name, args.out_name, args.chain_break_cutoff_A)


#/home/aljubetic/conda/envs/pyro/bin/python rechain.py /home/nbizjak/projects/11_04_2023_rigid_connections/tests/mpnn_stuff/01.pdb /home/nbizjak/projects/11_04_2023_rigid_connections/tests/mpnn_stuff/01.pdb 