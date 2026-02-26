import pickle
from Bio.PDB import PDBParser, PDBIO, Superimposer, PPBuilder
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio import SeqIO
from Bio.Seq import Seq
import json
import glob
import os
import numpy as np
import pandas as pd
import shutil
from pathlib import Path
import string
import homooligomer_rmsd
from Bio.Align import PairwiseAligner
import re
import yaml


def make_boltz_input_yaml(
    cfg, model_id, mpnn_sequence, output_dir, input_alignment_dir
):
    chain_ids = []
    sequences = []
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    # Split multi-chain sequences by colon
    print(
        f"splitting mpnn sequence {mpnn_sequence} by colon for Boltz yaml generation..."
    )
    split_chains = mpnn_sequence.split(":")
    for i, chain_seq in enumerate(split_chains):
        chain_id = f"{letters[i]}"
        chain_ids.append(chain_id)
        print(f"Chain ID: {chain_id}, Sequence: {chain_seq}")
        sequences.append(chain_seq)

    # Make boltz yaml
    data = dict(sequences=dict())
    data["sequences"] = []
    if cfg.use_a3m and input_alignment_dir is not None:
        for idxsequence, sequence in enumerate(sequences):
            chain_id = chain_ids[idxsequence]
            cleaned = "".join(c for c in sequence if c.isalpha())
            print(
                f"Chain ID: {chain_ids[idxsequence]}, Sequence: {sequence}, Cleaned Sequence: {cleaned}"
            )
            data["sequences"].append(
                {
                    "protein": {
                        "id": chain_id,
                        "sequence": cleaned,
                        "msa": f"{input_alignment_dir}/{model_id}_{chain_id}.a3m",
                    }
                }
            )
    else:
        for idxsequence, sequence in enumerate(sequences):
            cleaned = "".join(c for c in sequence if c.isalpha())
            print(
                f"Chain ID: {chain_ids[idxsequence]}, Sequence: {sequence}, Cleaned Sequence: {cleaned}"
            )
            chain_id = chain_ids[idxsequence]
            data["sequences"].append(
                {
                    "protein": {
                        "id": chain_id,
                        "sequence": cleaned,
                        "msa": "empty",
                    }
                }
            )

    with open(f"{output_dir}/{model_id}.yaml", "w") as outfile:
        yaml.dump(data, outfile, default_flow_style=False)
    return f"{output_dir}/{model_id}.yaml"


def get_rmsd_from_coords(native_coords, model_coords, rot, tran):
    model_coords_rotated = np.dot(model_coords, rot) + tran
    diff = native_coords - model_coords_rotated
    RMSD = np.sqrt(sum(sum(diff**2)) / native_coords.shape[0])
    return RMSD


def extract_chain_letter(filename):
    match = re.search(r"Chain_([A-Z])", filename)
    return match.group(1) if match else ""


def masked_positions(seq1, seq2, min_block=3):
    aligner = PairwiseAligner()
    aligner.mode = "local"

    alignment = aligner.align(seq1, seq2)[0]
    aligned1, aligned2 = alignment[0], alignment[1]

    masked = []
    positions = []
    seq1_index = 0

    for a, b in zip(aligned1, aligned2):
        if a != "-":
            masked.append(a == b)
            seq1_index += 1

    # Apply min_block filtering
    n = len(masked)
    i = 0
    while i < n:
        if not masked[i]:
            i += 1
            continue
        j = i
        while j < n and masked[j]:
            j += 1
        if (j - i) < min_block:
            for k in range(i, j):
                masked[k] = False
        i = j

    # Collect positions that remain True
    positions = [idx for idx, keep in enumerate(masked) if keep]

    return positions


def calculate_RMSD_linker_len(
    cfg, trb_path, af2_pdb, starting_pdb, rfdiff_pdb_path, symmetry, model_monomer
):
    # First calculate RMSD between input protein and AF2 generated protein
    # Second calcualte number of total generated AA by RFDIFF
    #   - if designing only in one location the number is equal linker length

    parser = PDBParser(PERMISSIVE=1)
    structure_af2 = parser.get_structure("af2", af2_pdb)

    skipRfDiff = cfg.get("skipRfDiff", False)
    designable_residues = cfg.get("designable_residues", None)
    # Skip if trb does not exist

    if skipRfDiff:
        print(
            "skipping rfdiffusion. RMSD_Sculpted is that of the designable residues. RMSD_Motif is that of the non-designable residues"
        )

        if symmetry != None or model_monomer:
            rmsd = homooligomer_rmsd.align_oligomers(
                starting_pdb, af2_pdb, save_aligned=False
            )
            return ([round(rmsd, 1), -1, -1, -1, -1], -1)
        else:  # There is a lot of code duplication here from when RFDiff is used, but I don't want to deal with the problem of making the skipped one go into the non-skip correctly
            # rmsd = homooligomer_rmsd.align_monomer(starting_pdb, af2_pdb, save_aligned=False)

            chainResidOffset, con_hal_pdb_idx_complete = getChainResidOffsets(
                cfg.pdb_path, designable_residues
            )

            selected_residue_data = [
                x + (chainResidOffset[chain] - 1)
                for chain, x in con_hal_pdb_idx_complete
            ]
            selected_residues_in_designed_chains = selected_residue_data  # This will be a problem if there are fixed chains. TODO:FIX

            all_af2_res = list(structure_af2.get_residues())
            all_af2_res_ca = [ind["CA"] for ind in all_af2_res]
            af2_all_fixed_res = [
                all_af2_res[ind]["CA"] for ind in selected_residue_data
            ]
            af2_sculpted_res = [
                ind["CA"] for ind in all_af2_res if ind["CA"] not in af2_all_fixed_res
            ]
            # af2_fixed_chain_res=[all_af2_res[ind]['CA'] for ind in selected_residues_in_fixed_chains]
            # af2_motif_res=[all_af2_res[ind]['CA'] for ind in selected_residues_in_designed_chains]
            structure_rfdiff = parser.get_structure("control", rfdiff_pdb_path)
            all_rfdiff_res = list(
                structure_rfdiff.get_residues()
            )  # obtain a list of all the residues in the structure, structure_control is object
            all_rfdiff_res_ca = [ind["CA"] for ind in all_rfdiff_res]

            superimposer = SVDSuperimposer()
            rfdiff_all_coords = np.array([a.coord for a in all_rfdiff_res_ca])
            af2_all_coords = np.array([a.coord for a in all_af2_res_ca])

            superimposer.set(rfdiff_all_coords, af2_all_coords)
            superimposer.run()
            rmsd = get_rmsd_from_coords(
                rfdiff_all_coords, af2_all_coords, superimposer.rot, superimposer.tran
            )

            rfdiff_all_fixed_res = [
                all_rfdiff_res[ind]["CA"] for ind in selected_residue_data
            ]  # retrieve the residue with the corresponding index from rfdiff_res
            rfdiff_sculpted_res = [
                ind["CA"]
                for ind in all_rfdiff_res
                if ind["CA"] not in rfdiff_all_fixed_res
            ]

            superimposer = SVDSuperimposer()
            rfdiff_all_fixed_coords = np.array([a.coord for a in rfdiff_all_fixed_res])
            af2_all_fixed_coords = np.array([a.coord for a in af2_all_fixed_res])
            superimposer.set(rfdiff_all_fixed_coords, af2_all_fixed_coords)
            superimposer.run()
            rmsd_all_fixed = get_rmsd_from_coords(
                rfdiff_all_fixed_coords,
                af2_all_fixed_coords,
                superimposer.rot,
                superimposer.tran,
            )
            # rfdiff_fixed_chain_res=[all_rfdiff_res[ind]['CA'] for ind in selected_residues_in_fixed_chains]
            # rfdiff_motif_res=[all_rfdiff_res[ind]['CA'] for ind in selected_residues_in_designed_chains]
            rfdiff_sculpted_coords = [a.coord for a in rfdiff_sculpted_res]
            af2_sculpted_coords = [a.coord for a in af2_sculpted_res]
            rfdiff_sculpted_coords = np.array(rfdiff_sculpted_coords)
            af2_sculpted_coords = np.array(af2_sculpted_coords)

            rmsd_sculpted = get_rmsd_from_coords(
                rfdiff_sculpted_coords,
                af2_sculpted_coords,
                superimposer.rot,
                superimposer.tran,
            )
            print(
                [
                    round(rmsd, 1),
                    round(rmsd_all_fixed, 1),
                    round(rmsd_sculpted, 1),
                    -1,
                    round(rmsd_all_fixed, 1),
                ]
            )
            return (
                [
                    round(rmsd, 1),
                    round(rmsd_all_fixed, 1),
                    round(rmsd_sculpted, 1),
                    -1,
                    round(rmsd_all_fixed, 1),
                ],
                -1,
            )

    with open(trb_path, "rb") as f:
        trb_dict = pickle.load(f)

    # Get different data from trb file depending on the fact if designing a monomer (one chain) or heteromer
    # complex_con_rex_idx present only if there are chains that are completely fixed.
    # Data structure: con_ref_idx0 = [0, 1, 3, 3, ...]
    #   Info: array of input pdb AA indices starting 0 (con_ref_pdb_idx), and where they are in the output pdb (con_hal_pdb_idx)
    #   In complex_con_hal_idx0 there is no chain info however RFDIFF changes pdb indeces to go from 1 to n (e.g. 1st AA in chain B has idx 34)

    # selected_residues_data will hold the information only of those residues that were selected from the reference structure to be used in the final design
    if "complex_con_ref_idx0" in trb_dict:
        selected_residues_data = trb_dict["complex_con_hal_idx0"]
        selected_residues_in_designed_chains = trb_dict["con_hal_idx0"]
        if (
            len(selected_residues_in_designed_chains) == 0
            and trb_dict["config"]["contigmap"]["provide_seq"] != None
        ):
            selected_residues_in_designed_chains = np.where(
                [
                    a != b
                    for a, b in zip(trb_dict["inpaint_seq"], trb_dict["inpaint_str"])
                ]
            )[
                0
            ]  # if this works...
            print(
                f"DEBUG: PARTIAL DIFUSSION KEEPING RESIDUES {selected_residues_in_designed_chains}"
            )
            # print(selected_residues_in_designed_chains)

        # selected_residues_in_fixed_chains=trb_dict['receptor_con_hal_idx0'] #This actually doesn't work and I think it's a bug in RFDiff
        selected_residues_in_fixed_chains = [
            res
            for res in selected_residues_data
            if res not in selected_residues_in_designed_chains
        ]
    else:
        partial_diffusion = cfg.get("partial_diffusion", False)
        if not partial_diffusion:
            selected_residues_data = trb_dict["con_hal_idx0"]
        else:  # When using partial diffusion, RFDiffusion doesn't put anything of the diffused chain on the
            # con_hal_pdb_idx (because nothing is technically fixed). This means that we need to recompile it
            # based on the inpaint_seq
            selected_residues_data = []

            chainResidOffset, con_hal_pdb_idx_complete = getChainResidOffsets(
                cfg.pdb_path, designable_residues
            )
            for id0, value in enumerate(trb_dict["inpaint_seq"]):
                if value == True:
                    abeceda = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                    for key in abeceda:
                        if key in chainResidOffset:
                            if id0 >= chainResidOffset[key]:
                                currentResidueChain = key
                    # residue_data_control_1.append((currentResidueChain,id0+chainResidOffset[currentResidueChain]))
                    selected_residues_data.append(id0)

        selected_residues_in_fixed_chains = []
        selected_residues_in_designed_chains = selected_residues_data

    all_af2_res = list(structure_af2.get_residues())
    all_af2_res_ca = [ind["CA"] for ind in all_af2_res]

    af2_all_fixed_res = [all_af2_res[ind]["CA"] for ind in selected_residues_data]
    af2_sculpted_res = [
        ind["CA"] for ind in all_af2_res if ind["CA"] not in af2_all_fixed_res
    ]
    af2_fixed_chain_res = [
        all_af2_res[ind]["CA"] for ind in selected_residues_in_fixed_chains
    ]
    af2_motif_res = [
        all_af2_res[ind]["CA"] for ind in selected_residues_in_designed_chains
    ]

    trb_help = list(trb_dict["inpaint_str"])
    linker_indeces = [
        boolean for boolean in trb_help if boolean == False
    ]  # calculate linker length here - convenient
    linker_length = len(linker_indeces)

    # io=PDBIO()
    # io.set_structure(structure_af2)
    # io.save("af2_pdb_2.pdb") #This is not necessary and might be slowing down everything a bit.
    rmsd = -1
    rmsd_all_fixed = (
        -1
    )  # If there's no starting structure, we cannot compare it. RMSD is undefined (-1)
    rmsd_sculpted = -1
    rmsd_fixed_chains = -1
    rmsd_motif = -1
    # if starting_pdb:
    structure_rfdiff = parser.get_structure("control", rfdiff_pdb_path)

    all_rfdiff_res = list(
        structure_rfdiff.get_residues()
    )  # obtain a list of all the residues in the structure, structure_control is object
    all_rfdiff_res_ca = [ind["CA"] for ind in all_rfdiff_res]

    rfdiff_all_fixed_res = [
        all_rfdiff_res[ind]["CA"] for ind in selected_residues_data
    ]  # retrieve the residue with the corresponding index from rfdiff_res
    rfdiff_sculpted_res = [
        ind["CA"] for ind in all_rfdiff_res if ind["CA"] not in rfdiff_all_fixed_res
    ]
    rfdiff_fixed_chain_res = [
        all_rfdiff_res[ind]["CA"] for ind in selected_residues_in_fixed_chains
    ]
    rfdiff_motif_res = [
        all_rfdiff_res[ind]["CA"] for ind in selected_residues_in_designed_chains
    ]

    if len(rfdiff_all_fixed_res) != len(af2_all_fixed_res):
        print(
            "Fixed and moving atom lists differ in size"
        )  # for now, this is when input pdb and output are different length
        print(rfdiff_all_fixed_res, af2_all_fixed_res)
        return (-1, -1)

    # Align all and get RMSD of all
    superimposer = SVDSuperimposer()
    rfdiff_all_coords = np.array([a.coord for a in all_rfdiff_res_ca])
    af2_all_coords = np.array([a.coord for a in all_af2_res_ca])

    superimposer.set(rfdiff_all_coords, af2_all_coords)
    superimposer.run()
    rmsd = get_rmsd_from_coords(
        rfdiff_all_coords, af2_all_coords, superimposer.rot, superimposer.tran
    )

    # Align all reference residues if there are no fixed chains. Otherwise, align only fixed chains.  (Very nice because fully fixed chains should be a stable reference)
    # then get rmsd_all_fixed and rmsd_sculpted (If any)
    # superimposer = SVDSuperimposer()
    rfdiff_all_fixed_coords = np.array([a.coord for a in rfdiff_all_fixed_res])
    af2_all_fixed_coords = np.array([a.coord for a in af2_all_fixed_res])

    rfdiff_fixed_chain_coords = np.array([a.coord for a in rfdiff_fixed_chain_res])
    rfdiff_motif_res_coords = np.array([a.coord for a in rfdiff_motif_res])

    af2_fixed_chain_coords = np.array([a.coord for a in af2_fixed_chain_res])
    af2_motif_res_coords = np.array([a.coord for a in af2_motif_res])

    if len(rfdiff_fixed_chain_coords) == 0:  # (there are no fixed chains)
        if len(rfdiff_all_fixed_coords) != 0:  # (There are fixed residues at all)
            superimposer.set(rfdiff_all_fixed_coords, af2_all_fixed_coords)
            superimposer.run()
            rmsd_all_fixed = get_rmsd_from_coords(
                rfdiff_all_fixed_coords,
                af2_all_fixed_coords,
                superimposer.rot,
                superimposer.tran,
            )
        else:
            rmsd_all_fixed = -1
    else:
        superimposer.set(rfdiff_fixed_chain_coords, af2_fixed_chain_coords)
        superimposer.run()
        rmsd_fixed_chains = get_rmsd_from_coords(
            rfdiff_fixed_chain_coords,
            af2_fixed_chain_coords,
            superimposer.rot,
            superimposer.tran,
        )
        rmsd_all_fixed = get_rmsd_from_coords(
            rfdiff_all_fixed_coords,
            af2_all_fixed_coords,
            superimposer.rot,
            superimposer.tran,
        )

    if True in trb_dict["inpaint_seq"]:  # There are non-redesigned residues
        rfdiff_sculpted_coords = [a.coord for a in rfdiff_sculpted_res]
        af2_sculpted_coords = [a.coord for a in af2_sculpted_res]

        rfdiff_sculpted_coords = np.array(rfdiff_sculpted_coords)
        af2_sculpted_coords = np.array(af2_sculpted_coords)

        # superimposer.set(rfdiff_all_coords, af2_all_coords)
        # superimposer.run()

        rmsd_sculpted = get_rmsd_from_coords(
            rfdiff_sculpted_coords,
            af2_sculpted_coords,
            superimposer.rot,
            superimposer.tran,
        )
        if len(rfdiff_motif_res_coords) != 0:
            rmsd_motif = get_rmsd_from_coords(
                rfdiff_motif_res_coords,
                af2_motif_res_coords,
                superimposer.rot,
                superimposer.tran,
            )

    # If we do symmetry, we align af2 model to rfdiffusion structure. Should we control that, or hardcode it?

    if symmetry != None:
        rmsd = homooligomer_rmsd.align_oligomers(
            rfdiff_pdb_path, af2_pdb, save_aligned=False
        )

    return (
        [
            round(rmsd, 1),
            round(rmsd_all_fixed, 1),
            round(rmsd_sculpted, 1),
            round(rmsd_fixed_chains, 1),
            round(rmsd_motif, 1),
        ],
        linker_length,
    )


def make_alignment_file_boltz(sequence_id, sequence, alignment_dir, output_dir):
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    name = sequence_id
    seq = sequence
    chain_ids = []
    sequences = []
    # Split multi-chain sequences by colon
    split_chains = seq.split(":")
    for i, chain_seq in enumerate(split_chains):
        chain_id = f"{letters[i]}"
        chain_ids.append(chain_id)
        sequences.append(chain_seq)

    chain_lengths = [len(s) for s in sequences]
    # Get all relevant files
    alignment_files = [f for f in os.listdir(alignment_dir) if f.endswith(".a3m")]

    # Sort by the chain letter
    alignment_files_sorted = sorted(alignment_files, key=extract_chain_letter)

    for idxsequence, sequence in enumerate(sequences):
        sequence = sequence.strip()
        chain_id = chain_ids[idxsequence]
        print(f"Processing chain {chain_id} for {name} with sequence {sequence}...")

        with open(f"{output_dir}/{name}_{chain_id}.a3m", "w+") as f:
            f.write(f">101\n")
            f.write(sequence + "\n")

            # Iterate in sorted order
            filename = None
            for fn in alignment_files_sorted:
                if f"Chain_{chain_id}" in fn:
                    filename = fn
                    break
            if filename is None:
                print(f"Could not find alignment file for chain {chain_id}")
                f.write(f">101\n")
                f.write(sequence + "\n")
            else:
                path = os.path.join(alignment_dir, filename)
                with open(path) as afile:

                    next(afile)
                    # read original sequence
                    reference_name = afile.readline()
                    reference_seq = afile.readline()
                    reference_seq = "".join(
                        [c for c in reference_seq if not c.islower()]
                    )
                    # Compute positions to keep
                    positions_to_keep = masked_positions(sequence, reference_seq)

                    while True:
                        name_line = afile.readline()
                        seq_line = afile.readline()
                        seq_line = "".join([c for c in seq_line if not c.islower()])
                        if not seq_line:
                            break

                        name_line = name_line.rstrip()
                        seq_line = seq_line.rstrip()

                        masked_seq = ["-"] * (len(sequence))  # NEEDS A -1
                        # Place masked residues at the correct positions in full_sequence
                        for j, pos in enumerate(positions_to_keep):
                            if j < len(seq_line):
                                masked_seq[pos] = seq_line[j]
                        f.write(name_line + "\n")
                        f.write("".join(masked_seq) + "\n")

        with open(f"{output_dir}/{name}_{chain_id}.a3m", "rb+") as f:
            f.seek(-1, 2)  # move to the last byte
            if f.read(1) == b"\n":
                f.seek(-1, 2)
                f.truncate()

        cleaned_lines = []

        # Read and clean
        with open(f"{output_dir}/{name}_{chain_id}.a3m", "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#") or line.startswith(">"):
                    cleaned_lines.append(line)
                    continue
                cleaned_seq = "".join([c.upper() if c.isalpha() else "-" for c in line])
                cleaned_lines.append(cleaned_seq)

        # Overwrite the same file
        with open(f"{output_dir}/{name}_{chain_id}.a3m", "w") as f:
            f.write("\n".join(cleaned_lines))


def make_alignment_file(cfg, trb_path, pdb_file, mpnn_seq, alignments_path, output):

    skipRfDiff = cfg.get("skipRfDiff", False)
    designable_residues = cfg.get("designable_residues", None)
    chainResidOffset, con_hal_pdb_idx_complete = getChainResidOffsets(
        pdb_file, designable_residues
    )

    if not skipRfDiff:
        with open(trb_path, "rb") as f:
            trb_dict = pickle.load(f)

        if "complex_con_ref_idx0" in trb_dict:
            # residue_data_control_0 = trb_dict['complex_con_ref_idx0']
            residue_data_af2_0 = trb_dict["complex_con_hal_idx0"]
            residue_data_control_1 = trb_dict["complex_con_ref_pdb_idx"]
            # residue_data_af2_1 = trb_dict['complex_con_hal_pdb_idx']
        else:

            partial_diffusion = cfg.get("partial_diffusion", False)
            if not partial_diffusion:
                residue_data_af2_0 = trb_dict["con_hal_idx0"]
                residue_data_control_1 = trb_dict["con_ref_pdb_idx"]
            else:  # When using partial diffusion, RFDiffusion doesn't put anything of the diffused chain on the
                # con_hal_pdb_idx (because nothing is technically fixed). This means that we need to recompile it
                # based on the inpaint_seq
                residue_data_control_1 = []
                residue_data_af2_0 = []

                for id0, value in enumerate(trb_dict["inpaint_seq"]):
                    if value == True:
                        abeceda = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                        for key in abeceda:
                            if key in chainResidOffset:
                                if id0 >= chainResidOffset[key]:
                                    currentResidueChain = key
                        residue_data_control_1.append(
                            (
                                currentResidueChain,
                                id0 + chainResidOffset[currentResidueChain],
                            )
                        )
                        residue_data_af2_0.append(id0)

            # residue_data_control_0 = trb_dict['con_ref_idx0']

            # residue_data_af2_1 = trb_dict['con_hal_pdb_idx']
    else:
        residue_data_af2_0 = [
            x + (chainResidOffset[chain] - 1) for chain, x in con_hal_pdb_idx_complete
        ]
        residue_data_control_1 = con_hal_pdb_idx_complete
        print(residue_data_af2_0)
        print(con_hal_pdb_idx_complete)

    if mpnn_seq[-1:] == "\n":
        mpnn_seq = mpnn_seq[:-1]
    mpnn_sequence_no_colons = mpnn_seq.replace(":", "")

    used_chains = list(set([i[0] for i in residue_data_control_1]))

    mpnn_sequences_list = mpnn_seq.split(":")
    sequences_limits = []

    for seq_num, sequence in enumerate(mpnn_sequences_list):
        seq_start = mpnn_sequence_no_colons.find(sequence)

        if seq_num != len(mpnn_sequences_list) - 1:
            next_sequence = mpnn_sequences_list[seq_num + 1]
            seq_end = mpnn_sequence_no_colons.find(next_sequence)
        else:
            seq_end = len(mpnn_sequence_no_colons)

        sequences_limits.append((seq_start, seq_end))

    with open(output, "w") as f:
        letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"  # used for naming chains

        # write the header line
        first_line = "#"
        for seq in mpnn_sequences_list:
            first_line += str(len(seq))
            if seq != mpnn_sequences_list[len(mpnn_sequences_list) - 1]:
                first_line += ","
        first_line += "\t"
        for seq in mpnn_sequences_list:
            first_line += "1"
            if seq != mpnn_sequences_list[len(mpnn_sequences_list) - 1]:
                first_line += ","
        f.write(first_line + "\n")

        # write the whole sequence once
        all_names = ""
        for seq_num, sequence in enumerate(mpnn_sequences_list):
            # all_names+=str(101+seq_num)
            all_names += letters[seq_num]
            if sequence != mpnn_sequences_list[len(mpnn_sequences_list) - 1]:
                all_names += "\t"

        f.write(">" + all_names + "\n")
        f.write(mpnn_sequence_no_colons + "\n")

        # write the sequences to be modelled:
        for seq_num, sequence in enumerate(mpnn_sequences_list):
            f.write(">" + letters[seq_num] + "\n")
            sequence_line = (
                "-" * sequences_limits[seq_num][0]
            )  # Add a gap for each position before the sequence
            sequence_line += sequence  # add sequence
            sequence_line += "-" * (
                (len(mpnn_sequence_no_colons) - sequences_limits[seq_num][1])
            )  # Add a gap for each position after the sequence.
            f.write(sequence_line + "\n")  # write padded sequence

        # now write the aligned sequences
        for chain in letters:
            if chain in used_chains:
                # LEt's get the correct file for this chain
                for file in os.listdir(alignments_path):
                    if "auth_" + chain in file or "Chain_" + chain in file:
                        alignment_file = file
                        print(
                            "Alignment file for chain "
                            + chain
                            + " is "
                            + alignment_file
                        )

                with open(
                    os.path.join(alignments_path, alignment_file), "r"
                ) as chain_alignment_file:
                    for line_id, line in enumerate(chain_alignment_file):
                        if (
                            line_id >= 3
                        ):  # skip first three lines, since they contain the original sequence.
                            if line[0] == ">":
                                f.write(line)
                            else:
                                table = str.maketrans(
                                    "", "", string.ascii_lowercase
                                )  # This deletes lowercase characters from the string
                                line_without_insertions = line.translate(table)

                                new_aligned_seq = "-" * (
                                    len(mpnn_sequence_no_colons)
                                )  # Make a gap sequence of the length of the sequence..
                                trb_chain = [
                                    x
                                    for x in residue_data_control_1
                                    if x[0][0] == chain
                                ]
                                first_residue_in_trb = trb_chain[0][1]
                                for id, pos in enumerate(residue_data_control_1):
                                    if (
                                        pos[0] == chain
                                    ):  # If position chain corresponds to the chain we're looking at

                                        position_to_copy = (
                                            residue_data_control_1[id][1] - 1
                                        )  # minus 1 because this is 1-indexed while the sequence is 0 indexed
                                        new_aligned_seq = (
                                            new_aligned_seq[: residue_data_af2_0[id]]
                                            + line_without_insertions[
                                                position_to_copy
                                                - first_residue_in_trb
                                                + 1
                                            ]
                                            + new_aligned_seq[
                                                residue_data_af2_0[id] + 1 :
                                            ]
                                        )

                                f.write(new_aligned_seq + "\n")

    # delete empty lines that are generated for weird reasons beyond my comprehension. This should be fixed and this section removed, but it doesn't really slow things that much.
    with open(output, "r+") as output_file:
        with open(output + "_tmp", "w") as temp_file:
            for line in output_file:
                if not line.isspace():
                    temp_file.write(line)

    os.remove(output)
    os.rename(output + "_tmp", output)
    # shutil.copyfile(output, output+"_backup") #this is for debug only, to see the file before it goes to AF2


def get_token_value(
    astr, token, regular_expression
):  # "(\d*\.\d+|\d+\.?\d*)" # (-?\d*\.\d+|-?\d+\.?\d*) to allow negative RMSD (-1 = undefined)
    """returns value next to token"""
    import re

    regexp = re.compile(f"{token}{regular_expression}")
    match = regexp.search(astr)
    # if match == None:
    #   match = "/"
    return match.group(1)


def merge_csv(output_dir, output_csv, scores_csv):
    # read csv files
    df1 = pd.read_csv(scores_csv)
    df2 = pd.read_csv(output_csv)

    # merge dataframes on 'model_path' column
    merged_df = pd.merge(df1, df2, on="model_path")
    # drop duplicate 'model_path' column (if it exists)
    merged_df = merged_df.loc[:, ~merged_df.columns.duplicated()]

    # save merged dataframe to csv file
    merged_df.to_csv(
        f'{os.path.join(output_dir, "final_output.csv")}',
        index=False,
        float_format="%.1f",
    )

    # Select best ones and copy to another csv. commented out for now
    # if (output_best):
    # best_df=merged_df[(merged_df["RMSD"] <= rmsd_threshold) | (merged_df["plddt"] >= plddt_threshold)] #select best based on thresholds
    # best_df.to_csv(f'{os.path.join(output_dir, "final_output_best.csv")}', index=False)

    # dir_best_pdbs = os.path.join(output_dir, "best_pdbs")
    # os.makedirs(dir_best_pdbs, exist_ok=True) # directory is created even if some or all of the intermediate directories in the path do not exist

    # for file in best_df["model_path"]:
    #    shutil.copy(file,os.path.join(output_dir, "best_pdbs"))


def rename_pdb_create_csv_colabfold(
    cfg,
    output_dir,
    rfdiff_out_dir,
    trb_num,
    model_i,
    control_structure_path,
    symmetry=None,
    model_monomer=False,
):

    # Preparing paths to acces correct files
    model_i = os.path.join(model_i, "")  # add / to path to access json files within

    # dir_renamed_pdb = os.path.join(os.path.dirname(output_dir), "final_pdbs") #Why is this done to the parent folder? It's annoying if running multiple jobs on the same folder
    dir_renamed_pdb = os.path.join(output_dir, "final_pdbs")
    os.makedirs(
        dir_renamed_pdb, exist_ok=True
    )  # directory is created even if some or all of the intermediate directories in the path do not exist

    trb_file = os.path.join(
        rfdiff_out_dir, f"_{trb_num}.trb"
    )  # name of corresponding trb file
    skipRfDiff = cfg.get("skipRfDiff", False)
    designable_residues = cfg.get("designable_residues", None)
    partial_diffusion = cfg.get("partial_diffusion", False)

    if not skipRfDiff:
        with open(trb_file, "rb") as f:
            trb_dict = pickle.load(f)

        if not partial_diffusion:
            if "complex_con_ref_idx0" in trb_dict:
                residue_data_af2 = trb_dict["complex_con_hal_idx0"]
            else:
                residue_data_af2 = trb_dict["con_hal_idx0"]
        else:  # When using partial diffusion, RFDiffusion doesn't put anything of the diffused chain on the
            # con_hal_pdb_idx (because nothing is technically fixed). This means that we need to recompile it
            # based on the inpaint_seq
            residue_data_af2 = []
            abeceda = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            chainResidOffset, con_hal_pdb_idx_complete = getChainResidOffsets(
                control_structure_path, designable_residues
            )
            for id0, value in enumerate(trb_dict["inpaint_seq"]):

                if value == True:

                    for key in abeceda:
                        if key in chainResidOffset:

                            if id0 >= chainResidOffset[key]:
                                currentResidueChain = key
                    residue_data_af2.append(id0)
                    # residue_data_af2.append((currentResidueChain,id0+chainResidOffset[currentResidueChain]))

    else:
        chainResidOffset, con_hal_pdb_idx_complete = getChainResidOffsets(
            control_structure_path, designable_residues
        )
        residue_data_af2 = [
            x + (chainResidOffset[chain] - 1) for chain, x in con_hal_pdb_idx_complete
        ]

    json_files = glob.glob(os.path.join(model_i, "T*000.json"))
    rfdiff_pdb_path = os.path.join(rfdiff_out_dir, f"_{trb_num}.pdb")

    for (
        json_file
    ) in (
        json_files
    ):  # in af2 model_i directory for T_...json file in [all T_...json files]
        # This is done for each model_i directory therefore for each rfdiff pdb
        # There are 5 T_...jsons per 1 mpnn seq of the rfdiff model

        with open(json_file, "r") as f:
            params = json.load(f)

        # Handle filenames correctly to get the T_...pdb that corresponds to the T_...json
        json_filename = os.path.basename(json_file)
        json_dirname = os.path.dirname(json_file)
        json_newname = os.path.join(
            json_dirname, json_filename.replace("scores", "unrelaxed")
        )
        model_pdb_file = (
            os.path.splitext(json_newname)[0] + ".pdb"
        )  # T_0.1__sample_2__score_0.5830__global_score_0.8339__seq_recovery_02000_unrelaxed_rank_004_alphafold2_multimer_v3_model_3_seed_000.pdb

        # Extract relevant data. Files used: json file of specific af2 model, specific af2 pdb,  trb file of rfdiff model (1 for all AF2 models from same rfdiff pdb)
        plddt_list = params["plddt"]
        plddt = int(np.mean(plddt_list))

        try:
            plddt_sculpted_list = [
                plddt_list[i]
                for i in range(0, len(plddt_list))
                if i not in residue_data_af2
            ]

            plddt_sculpted = int(np.mean(plddt_sculpted_list))
        except NameError:
            plddt_sculpted = -1

        rmsd_list, linker_length = calculate_RMSD_linker_len(
            cfg,
            trb_file,
            model_pdb_file,
            control_structure_path,
            rfdiff_pdb_path,
            symmetry,
            model_monomer,
        )
        pae = round((np.mean(params["pae"])), 2)

        # if we are doing symmetry or monomer modelling we also want to add monomer rmsd to the output
        if symmetry:
            monomers_dirname = os.path.join(model_i, "monomers")
            monomer_pdb_file = os.path.join(
                monomers_dirname, "monomer_" + os.path.basename(model_pdb_file)
            )
            monomer_rmsd = homooligomer_rmsd.align_monomer(
                rfdiff_pdb_path, monomer_pdb_file, save_aligned=False
            )
            monomer_params_json = os.path.join(
                monomers_dirname, "monomer_" + os.path.basename(json_file)
            )
            with open(monomer_params_json, "r") as f:
                monomer_params = json.load(f)

            monomer_plddt_list = monomer_params["plddt"]
            monomer_plddt = int(np.mean(monomer_plddt_list))

        if model_monomer:
            monomers_dirname = os.path.join(model_i, "monomers")
            monomer_pdb_file = os.path.join(
                monomers_dirname, "monomer_" + os.path.basename(model_pdb_file)
            )

            parser = PDBParser(PERMISSIVE=1)

            structure_target = parser.get_structure("target", rfdiff_pdb_path)
            structure_mobile = parser.get_structure("mobile", monomer_pdb_file)

            target_chain = list(structure_target.get_chains())[0]
            mobile_chain_res = list(structure_mobile.get_residues())
            mobile_chain_res = [ind["CA"] for ind in mobile_chain_res]
            list_rmsd_chains = []

            target_chain_res = list(target_chain.get_residues())
            target_chain_res = [ind["CA"] for ind in target_chain_res]

            superimposer = Superimposer()
            superimposer.set_atoms(target_chain_res, mobile_chain_res)
            superimposer.apply(structure_mobile.get_atoms())
            list_rmsd_chains.append(superimposer.rms)
            monomer_rmsd = np.min(list_rmsd_chains)
            monomer_params_json = os.path.join(
                monomers_dirname, "monomer_" + os.path.basename(json_file)
            )
            with open(monomer_params_json, "r") as f:
                monomer_params = json.load(f)

            monomer_plddt_list = monomer_params["plddt"]
            monomer_plddt = int(np.mean(monomer_plddt_list))

        # tracebility
        output_num = os.path.basename(output_dir)
        af2_model = get_token_value(
            json_filename, "_model_", "(\\d*\\.\\d+|\\d+\\.?\\d*)"
        )
        mpnn_sample = get_token_value(
            json_filename, "_sample_", "(\\d*\\.\\d+|\\d+\\.?\\d*)"
        )
        task_id = os.environ.get("SLURM_ARRAY_TASK_ID", 1)

        # Create a new name an copy te af2 model under that name into the output directory
        new_pdb_file = f"{task_id}.{trb_num}.{mpnn_sample}.{af2_model}__link_{linker_length}__plddt_{plddt}__plddt_sculpted_{plddt_sculpted}__rmsd_{rmsd_list[0]:.1f}__rmsd_sculpted_{rmsd_list[2]:.1f}__rmsd_fixedchains_{rmsd_list[3]:.1f}__rmsd_motif_{rmsd_list[4]:.1f}__pae_{pae}__out_{output_num}_.pdb"
        # out -> 00 -> number of task
        # rf -> 01 -> number of corresponding rf difff model
        # af_model -> 4 -> number of the af model (1-5), can be set using --model_order flag
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

        print("new_pdb_file", new_pdb_file)
        dictionary = {
            "id": f"{task_id}.{trb_num}.{mpnn_sample}.{af2_model}",
            "link_lenght": (linker_length),
            "plddt": (plddt),
            "plddt_sculpted": (plddt_sculpted),
            "RMSD": f"{rmsd_list[0]:.1f}",
            #'Rmsd_all_fixed': get_token_value(new_pdb_file, '__rmsd_all_fixed_', "(-?\\d*\\.\\d+|-?\\d+\\.?\\d*)"),
            "RMSD_sculpted": f"{rmsd_list[2]:.1f}",
            "RMSD_fixed_chains": f"{rmsd_list[3]:.1f}",
            "RMSD_motif": f"{rmsd_list[4]:.1f}",
            "pae": (pae),
            "model_path": new_pdb_path,
            "sequence": seq[1:],
            "af2_json": json_file,
            "af2_pdb": model_pdb_file,
            "path_rfdiff": rfdiff_pdb_path,
        }  # MODEL PATH for scoring_rg_... #jsonfilename for traceability

        if symmetry or model_monomer:
            dictionary["monomer_rmsd"] = monomer_rmsd
            dictionary["monomer_plddt"] = monomer_plddt

        df = pd.json_normalize(dictionary)
        path_csv = os.path.join(output_dir, "output.csv")
        df.to_csv(
            path_csv,
            mode="a",
            header=not os.path.exists(path_csv),
            index=False,
            float_format="%.1f",
        )


def rename_pdb_create_csv_boltz(
    cfg,
    output_dir,
    rfdiff_out_dir,
    trb_num,
    model_i,
    control_structure_path,
    symmetry=None,
    model_monomer=False,
):

    # Preparing paths to acces correct files
    model_i = os.path.join(model_i, "")  # add / to path to access json files within

    # dir_renamed_pdb = os.path.join(os.path.dirname(output_dir), "final_pdbs") #Why is this done to the parent folder? It's annoying if running multiple jobs on the same folder
    dir_renamed_pdb = os.path.join(output_dir, "final_pdbs")
    os.makedirs(
        dir_renamed_pdb, exist_ok=True
    )  # directory is created even if some or all of the intermediate directories in the path do not exist

    trb_file = os.path.join(
        rfdiff_out_dir, f"_{trb_num}.trb"
    )  # name of corresponding trb file
    skipRfDiff = cfg.get("skipRfDiff", False)
    designable_residues = cfg.get("designable_residues", None)
    partial_diffusion = cfg.get("partial_diffusion", False)

    if not skipRfDiff:
        with open(trb_file, "rb") as f:
            trb_dict = pickle.load(f)

        if not partial_diffusion:
            if "complex_con_ref_idx0" in trb_dict:
                residue_data_af2 = trb_dict["complex_con_hal_idx0"]
            else:
                residue_data_af2 = trb_dict["con_hal_idx0"]
        else:  # When using partial diffusion, RFDiffusion doesn't put anything of the diffused chain on the
            # con_hal_pdb_idx (because nothing is technically fixed). This means that we need to recompile it
            # based on the inpaint_seq
            residue_data_af2 = []
            abeceda = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            chainResidOffset, con_hal_pdb_idx_complete = getChainResidOffsets(
                control_structure_path, designable_residues
            )
            for id0, value in enumerate(trb_dict["inpaint_seq"]):

                if value == True:

                    for key in abeceda:
                        if key in chainResidOffset:

                            if id0 >= chainResidOffset[key]:
                                currentResidueChain = key
                    residue_data_af2.append(id0)
                    # residue_data_af2.append((currentResidueChain,id0+chainResidOffset[currentResidueChain]))

    else:
        chainResidOffset, con_hal_pdb_idx_complete = getChainResidOffsets(
            control_structure_path, designable_residues
        )
        residue_data_af2 = [
            x + (chainResidOffset[chain] - 1) for chain, x in con_hal_pdb_idx_complete
        ]

    individual_directories = glob.glob(
        os.path.join(model_i, "boltz_results_yaml_inputs", "predictions", "*")
    )
    rfdiff_pdb_path = os.path.join(rfdiff_out_dir, f"_{trb_num}.pdb")

    for (
        directory
    ) in individual_directories:  # for each directory in the predictions folder
        dir_path = Path(directory)
        model_name = dir_path.parent.name if not dir_path.is_dir() else dir_path.name

        json_file = glob.glob(os.path.join(directory, "*.json"))[0]
        pae_file = glob.glob(os.path.join(directory, "pae*.npz"))[0]
        pde_file = glob.glob(os.path.join(directory, "pde*.npz"))[0]
        plddt_file = glob.glob(os.path.join(directory, "plddt*.npz"))[0]

        with open(json_file, "r") as f:
            params = json.load(f)
        plddt_list = np.load(plddt_file)["plddt"].tolist()
        plddt = int(np.mean(plddt_list) * 100)
        pae_list = np.load(pae_file)["pae"].tolist()
        pae = np.mean(pae_list)
        pde_list = np.load(pde_file)["pde"].tolist()
        pde = np.mean(pde_list)

        model_pdb_file = glob.glob(os.path.join(directory, "*.pdb"))[0]
        # print(f"DEBUG:residue_data_af2 {residue_data_af2}")
        try:
            plddt_sculpted_list = [
                plddt_list[i]
                for i in range(0, len(plddt_list))
                if i not in residue_data_af2
            ]

            plddt_sculpted = int(np.mean(plddt_sculpted_list) * 100)
        except NameError:
            plddt_sculpted = -1

        rmsd_list, linker_length = calculate_RMSD_linker_len(
            cfg,
            trb_file,
            model_pdb_file,
            control_structure_path,
            rfdiff_pdb_path,
            symmetry,
            model_monomer,
        )

        # if we are doing symmetry or monomer modelling we also want to add monomer rmsd to the output
        if symmetry:
            monomers_dirname = os.path.join(model_i, "monomers")
            monomer_pdb_file = os.path.join(
                monomers_dirname,
                "boltz_results_yaml_inputs",
                "predictions",
                "monomer_" + model_name,
                "monomer_" + os.path.basename(model_pdb_file),
            )
            monomer_rmsd = homooligomer_rmsd.align_monomer(
                rfdiff_pdb_path, monomer_pdb_file, save_aligned=False
            )
            monomer_json_file = glob.glob(os.path.join(directory, "*.json"))[0]
            monomer_pae_file = glob.glob(os.path.join(directory, "pae*.npz"))[0]
            monomer_pde_file = glob.glob(os.path.join(directory, "pde*.npz"))[0]
            monomer_plddt_file = glob.glob(os.path.join(directory, "plddt*.npz"))[0]
            with open(monomer_json_file, "r") as f:
                monomer_params = json.load(f)
            monomer_plddt_list = np.load(monomer_plddt_file)["plddt"].tolist()
            monomer_plddt = int(np.mean(monomer_plddt_list) * 100)
            monomer_pae_list = np.load(monomer_pae_file)["pae"].tolist()
            monomer_pae = np.mean(monomer_pae_list)
            monomer_pde_list = np.load(monomer_pde_file)["pde"].tolist()
            monomer_pde = np.mean(monomer_pde_list)

        if model_monomer:
            monomers_dirname = os.path.join(model_i, "monomers")
            monomer_pdb_file = os.path.join(
                monomers_dirname,
                "boltz_results_yaml_inputs",
                "predictions",
                "monomer_" + model_name,
                "monomer_" + os.path.basename(model_pdb_file),
            )

            parser = PDBParser(PERMISSIVE=1)

            structure_target = parser.get_structure("target", rfdiff_pdb_path)
            structure_mobile = parser.get_structure("mobile", monomer_pdb_file)

            target_chain = list(structure_target.get_chains())[0]
            mobile_chain_res = list(structure_mobile.get_residues())
            mobile_chain_res = [ind["CA"] for ind in mobile_chain_res]
            list_rmsd_chains = []

            target_chain_res = list(target_chain.get_residues())
            target_chain_res = [ind["CA"] for ind in target_chain_res]

            superimposer = Superimposer()
            superimposer.set_atoms(target_chain_res, mobile_chain_res)
            superimposer.apply(structure_mobile.get_atoms())
            list_rmsd_chains.append(superimposer.rms)
            monomer_rmsd = np.min(list_rmsd_chains)
            monomer_json_file = glob.glob(os.path.join(directory, "*.json"))[0]
            monomer_pae_file = glob.glob(os.path.join(directory, "pae*.npz"))[0]
            monomer_pde_file = glob.glob(os.path.join(directory, "pde*.npz"))[0]
            monomer_plddt_file = glob.glob(os.path.join(directory, "plddt*.npz"))[0]
            with open(monomer_json_file, "r") as f:
                monomer_params = json.load(f)
            monomer_plddt_list = np.load(monomer_plddt_file)["plddt"].tolist()
            monomer_plddt = int(np.mean(monomer_plddt_list) * 100)
            monomer_pae_list = np.load(monomer_pae_file)["pae"].tolist()
            monomer_pae = np.mean(monomer_pae_list)
            monomer_pde_list = np.load(monomer_pde_file)["pde"].tolist()
            monomer_pde = np.mean(monomer_pde_list)

        # tracebility
        output_num = os.path.basename(output_dir)
        af2_model = "1"  # Boltz only runs one model, so we can just set this to 1. If we wanted to run multiple AF2 models per RFDiffusion model, we would need to change the Boltz code to output the model number in the json filename and then extract it here like we do for the regular ColabFold runs.
        mpnn_sample = get_token_value(
            model_name, "_sample_", "(\\d*\\.\\d+|\\d+\\.?\\d*)"
        )
        task_id = os.environ.get("SLURM_ARRAY_TASK_ID", 1)

        # Create a new name an copy te af2 model under that name into the output directory
        new_pdb_file = f"{task_id}.{trb_num}.{mpnn_sample}.{af2_model}__link_{linker_length}__plddt_{plddt}__plddt_sculpted_{plddt_sculpted}__rmsd_{rmsd_list[0]:.1f}__rmsd_sculpted_{rmsd_list[2]:.1f}__rmsd_fixedchains_{rmsd_list[3]:.1f}__rmsd_motif_{rmsd_list[4]:.1f}__pae_{pae}__out_{output_num}_.pdb"
        # out -> 00 -> number of task
        # rf -> 01 -> number of corresponding rf difff model
        # af_model -> 4 -> number of the af model (1-5), can be set using --model_order flag
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

        print("new_pdb_file", new_pdb_file)
        dictionary = {
            "id": f"{task_id}.{trb_num}.{mpnn_sample}.{af2_model}",
            "link_lenght": (linker_length),
            "plddt": (plddt),
            "plddt_sculpted": (plddt_sculpted),
            "RMSD": f"{rmsd_list[0]:.1f}",
            #'Rmsd_all_fixed': get_token_value(new_pdb_file, '__rmsd_all_fixed_', "(-?\\d*\\.\\d+|-?\\d+\\.?\\d*)"),
            "RMSD_sculpted": f"{rmsd_list[2]:.1f}",
            "RMSD_fixed_chains": f"{rmsd_list[3]:.1f}",
            "RMSD_motif": f"{rmsd_list[4]:.1f}",
            "pae": (pae),
            "model_path": new_pdb_path,
            "sequence": seq[1:],
            "af2_json": json_file,
            "af2_pdb": model_pdb_file,
            "path_rfdiff": rfdiff_pdb_path,
        }  # MODEL PATH for scoring_rg_... #jsonfilename for traceability
        dictionary.update(params)
        if symmetry or model_monomer:
            dictionary["monomer_rmsd"] = monomer_rmsd
            dictionary["monomer_plddt"] = monomer_plddt

        df = pd.json_normalize(dictionary)
        path_csv = os.path.join(output_dir, "output.csv")
        df.to_csv(
            path_csv,
            mode="a",
            header=not os.path.exists(path_csv),
            index=False,
            float_format="%.1f",
        )


# vsaka tabelca za svoj task
# funkcija na koncu, ki vse združi
# sestavljanje pathov je ok, dodaj nov column s pathom do rfdif


def create_dataframe(path_to_files, output_dir):  # path = r'content/*partial.pdb'
    # takes path  to renamed pdbs
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
        dictionary = {
            "link_lenght": get_token_value(
                file_name, "link_", "(\\d*\\.\\d+|\\d+\\.?\\d*)"
            ),
            "plddt": get_token_value(
                file_name, "__plddt_", "(\\d*\\.\\d+|\\d+\\.?\\d*)"
            ),
            "RMSD": get_token_value(
                file_name, "__rmsd_", "(-?\\d*\\.\\d+|-?\\d+\\.?\\d*)"
            ),
            #'Rmsd_all_fixed': get_token_value(file_name, '__rmsd_all_fixed_', "(-?\\d*\\.\\d+|-?\\d+\\.?\\d*)"),
            "RMSD_sculpted": get_token_value(
                file_name, "__rmsd_sculpted_", "(-?\\d*\\.\\d+|-?\\d+\\.?\\d*)"
            ),
            "pae": get_token_value(file_name, "__pae_", "(\\d*\\.\\d+|\\d+\\.?\\d*)"),
            "model_path": file_name,
            "sequence": seq,
            "rfdiff model": get_token_value(
                file_name, "__rf_", "(\\d*\\.\\d+|\\d+\\.?\\d*)"
            ),
        }  # MODEL PATH for scoring_rg_... #jsonfilename for traceability

        list_of_dicts.append(dictionary)

    # columns = ['link_length', 'plddt', 'loop_plddt', 'RMSD', 'model_path', 'sequence', 'score_traceb']
    df = pd.DataFrame(list_of_dicts)
    path_csv = os.path.join(os.path.dirname(output_dir), "output.csv")
    df.to_csv(
        path_csv,
        mode="a",
        header=not os.path.exists(path_csv),
        index=False,
        float_format="%.1f",
    )


class NumpyInt64Encoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.int64):
            return int(obj)
        return super(NumpyInt64Encoder, self).default(obj)


def getChainResidOffsets(pdb_file, designable_residues):
    chainResidOffset = {}
    con_hal_idx = []

    parser = PDBParser()
    chainsInPdb = parser.get_structure("protein", pdb_file).get_chains()
    for let_chain in chainsInPdb:
        chainLetter = let_chain.get_id()
        let_resids = let_chain.get_residues()
        startingNo = int(next(let_resids).get_id()[1])
        chainResidOffset.setdefault(chainLetter, startingNo - 1)
        if designable_residues:
            # We do this here to avoid PDBparser overhead
            print("There is designable_residues. [inside getChainResidOffsets]")
            for r in let_chain.get_residues():
                # aa has id " ". Heteroatoms have id "W" for water etc. – we don't want them in the PDB (they count as residues and ProteinMPNN throws an out-of-range error)
                if (
                    r.get_id()[0].strip() == ""
                    and f"{chainLetter}{r.get_id()[1]}" not in designable_residues
                ):
                    con_hal_idx.append((chainLetter, r.get_id()[1]))
    return chainResidOffset, con_hal_idx


def process_pdb_files(pdb_path: str, out_path: str, cfg, trb_paths=None, cycle=0):
    skipRfDiff = cfg.get("skipRfDiff", False)
    designable_residues = cfg.get("designable_residues", None)

    fixpos = {}
    pdb_files = Path(pdb_path).glob("*.pdb")

    contig = cfg.contig
    abeceda = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    for pdb_file in pdb_files:

        pdb_basename = pdb_file.stem

        ## We get RfDiff model number
        rf_model_num = get_token_value(
            os.path.basename(pdb_file), "_" if cycle == 0 else "rf_", "(\\d+)"
        )  # get 0 from _0.fa using reg exp #get 0 from rf_0__model_1__cycle_2__itr_0__.pdb
        # (We could also just get the first _# and exit regex early -- RfDiff model number is always the first element of the filename)
        print(f"RfDiff model number: {rf_model_num}")

        fixed_res = {}

        # We need to renumber fixed resids: each chain should start with 1
        chainResidOffset, con_hal_pdb_idx_complete = getChainResidOffsets(
            pdb_file, designable_residues
        )
        print(f"ChainResidOffset: {chainResidOffset}")

        with open(
            f"{pdb_path}/../chainResidOffset_{rf_model_num}.json",
            "w" if cycle == 0 else "r",
        ) as f:
            if cycle == 0:
                json.dump(chainResidOffset, f)
            else:
                chainResidOffset = json.load(f)
        print(f"Dumped or read chainResidOffset_{rf_model_num}.json")
        print(f"ChainResidOffset: {chainResidOffset}")

        if not skipRfDiff:
            trb_file = pdb_file.with_suffix(".trb")

            if not trb_file.exists():
                rf_model_num = get_token_value(
                    os.path.basename(pdb_file), "rf_", "(\\d+)"
                )
                trb_file = os.path.join(
                    os.path.dirname(pdb_file), f"_{rf_model_num}.trb"
                )
                trb_file = trb_file.replace("2_1_cycle_directory", "1_rfdiff")
                print(
                    f"TRB file not found for {pdb_basename}. CAUTION, using composed path"
                )

            print(pdb_file, trb_file)

            with open(trb_file, "rb") as f:
                trb_data = pickle.load(f)

            contig = trb_data["config"]["contigmap"]["contigs"][0]

            if "complex_con_hal_pdb_idx" in trb_data:
                con_hal_idx = trb_data.get(
                    "complex_con_hal_pdb_idx", []
                )  # con_hal_pdb_idx #complex_con_hal_pdb_idx
            else:
                partial_diffusion = cfg.get("partial_diffusion", False)
                if not partial_diffusion:
                    con_hal_idx = trb_data.get("con_hal_pdb_idx", [])
                else:  # When using partial diffusion, RFDiffusion doesn't put anything of the diffused chain on the
                    # con_hal_pdb_idx (because nothing is technically fixed). This means that we need to recompile it
                    # based on the inpaint_seq
                    con_hal_idx = []
                    for id0, value in enumerate(trb_data["inpaint_seq"]):
                        if value == True:
                            for key in abeceda:
                                if key in chainResidOffset:

                                    if id0 >= chainResidOffset[key]:
                                        currentResidueChain = key
                            con_hal_idx.append(
                                (
                                    currentResidueChain,
                                    id0 + chainResidOffset[currentResidueChain],
                                )
                            )

            # Process con_hal_idx to extract chain ids and indices
        else:
            con_hal_idx = con_hal_pdb_idx_complete

        # fixed_res should never be empty, otherwise ProteinMPNN will throw a KeyError fixed_position_dict[b['name']][letter]
        # We need to set blank fixed_res for each generated chain (based on contig).

        if cfg.inference.symmetry != None:  # if we find symmetry:
            breaks = int(cfg.inference.symmetry[1:])
        else:
            breaks = contig.count("/0 ") + 1

        fixed_res = dict(zip(abeceda, [[] for _ in range(breaks)]))
        print(f"Fixed res (according to contig chain breaks): {fixed_res}")

        # This is only good if multiple chains due to symmetry: all of them are equal; ProteinMPNN expects fixed_res as 1-based, resetting for each chain.

        for chain, idx in con_hal_idx:
            # If there are multiple chains, reset the auto_incrementing numbers to 1 for each chain (subtract offset)
            if not skipRfDiff:
                if trb_data["inpaint_seq"][
                    idx - 1
                ]:  # skip residues with FALSE in the inpaint_seq array
                    fixed_res.setdefault(chain, list()).append(
                        idx - chainResidOffset[chain]
                    )
            else:
                fixed_res.setdefault(chain, list()).append(
                    idx - chainResidOffset[chain]
                )
            # RfDiff outputs multiple chains if contig has /0 (chain break)

        print(f"Fixed res: ${fixed_res}")

        fixpos[pdb_basename] = fixed_res

    # print("_________trb data____", trb_data)

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
        seqs = Seq(sequence)
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
        # sequences = list(SeqIO.parse(file_handle, "fasta")) # better, for loop ni potreben
    return sequences


def change_sequence_in_fasta(pdb_file, mpnn_fasta):
    """
        #function adds : for chainbrakes and chain sequences of all other chains
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
    print(os.path.exists(pdb_file))
    i = 0
    seq_dict = {}
    for record in SeqIO.parse(mpnn_fasta, "fasta"):
        if i == 0:
            # Skip if record.description doesn't contain sample=. First seq is actually input to mpnn. However, after restart, we should not remove it again.
            if "sample=" not in record.description:
                print(
                    f"Skipping {record.description} for it does not contain sample=, it is input to and not output of mpnn."
                )
                i += 1
                continue
        i += 1
        print(record.seq, record.description)
        seq_dict[record.seq] = record.description

    print(seq_dict)
    sequences = []
    for record in SeqIO.parse(mpnn_fasta, "fasta"):
        if record.description in seq_dict.values():
            newseq = record.seq.replace("/", ":")
            record.seq = newseq
            sequences.append(record)
    print(sequences)
    SeqIO.write(
        sequences, mpnn_fasta, "fasta-2line"
    )  # This needs to be fasta-2line for MSA code to work


def match_linker_length(trb_path):
    import re

    with open(trb_path, "rb") as f:
        trb_dict = pickle.load(f)

    input_string = " ".join(
        trb_dict["sampled_mask"]
    )  #' '.join(['A1-30/6-6/C1-30/0', 'D1-30/0', 'B1-30/0'])

    pattern = r"(?<=/)(\\d+)(?=-\\d+/)"
    """
        (?<=/): Positive lookbehind assertion -> sub-pattern is preceded by a '/'.
        (\\d+): Capture group that matches one or more digits.
        (?=-\\d+/): Positive lookahead assertion -> sub-pattern is followed by a '-' and one or more digits, then a '/'.
    """
    match = re.search(pattern, input_string)
    if match:
        link = match.group(1)
        return link  # Output: 6
    else:
        print("No match found in linker length function")


def find_correct_trb(model_i):
    json_files = glob.glob(os.path.join(model_i, "_*000.json"))

    print(json_files)
    for json_file in json_files:
        json_filename = os.path.basename(json_file)
        number_trb = get_token_value(
            json_filename, "_", "(\\d+)(?=_)"
        )  # stop at first digit -> (?=_)
        trb_file = os.path.join(f"_{number_trb}.trb")
    return trb_file
