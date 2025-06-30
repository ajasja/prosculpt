# plugins/clash_check.py
from pyrosetta import *
from pyrosetta import Pose, pose_from_file
from pyrosetta.rosetta.protocols.toolbox.pose_manipulation import superimpose_pose_on_subset_CA
from pyrosetta.rosetta.utility import vector1_unsigned_long
from pyrosetta.rosetta.numeric import xyzVector_double_t
from scipy.spatial import cKDTree
import numpy as np
import os
import time
from pathlib import Path
from plugin_utils import ensure_pyrosetta_initialized
import logging

logger = logging.getLogger("prosculpt.plugins")

def make_vector1_range(start, end):
    vec = vector1_unsigned_long()
    for i in range(start, end + 1):
        vec.append(i)
    return vec

def copy_pose(pose):
    new_pose = Pose()
    new_pose.assign(pose)
    return new_pose

def get_chain_residues(pose, chain_id):
    return [i for i in range(1, pose.total_residue() + 1)
            if pose.pdb_info().chain(i) == chain_id]

def get_chain_residue_indices(pose, chain_index, residues_per_chain):
    start = 1 + chain_index * residues_per_chain
    end = start + residues_per_chain - 1
    return range(start, end + 1)

def is_valid_xyz(vec):
    return isinstance(vec, xyzVector_double_t) and hasattr(vec, 'x') and hasattr(vec, 'y') and hasattr(vec, 'z')

def assemble_c8_design(base_pdb, design_pdb, residues_to_align):
    base = pose_from_file(base_pdb)
    design = pose_from_file(design_pdb)

    design_chain_A_residues = get_chain_residues(design, 'A')
    if len(design_chain_A_residues) < residues_to_align:
        raise ValueError("Design chain A is too short for alignment.")

    design_align_residues = design_chain_A_residues[-residues_to_align:]
    combined = None

    for i in range(8):
        base_start = 1 + i * 301
        base_align_vec = make_vector1_range(base_start + 301 - residues_to_align, base_start + 300)
        aligned_copy = copy_pose(design)
        offset = -1 * ((base_start + 301 - residues_to_align) - design_align_residues[0])

        superimpose_pose_on_subset_CA(aligned_copy, base, base_align_vec, offset=offset)

        if combined is None:
            combined = copy_pose(aligned_copy)
        else:
            combined.append_pose_by_jump(aligned_copy, 1)
    
    # Dump the combined pose to a PDB file
    output_dir = Path("/home/d12-studenti/Kiyan/prosculpt/test_plugins/closest_pairs")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"assembled_{Path(design_pdb).stem}.pdb"
    combined.dump_pdb(str(output_path))
    logger.info(f"Dumped assembled pose to: {output_path}")
            
    return combined

def has_chain_pair_clashes(pose, chain_i, chain_j, residues_per_chain, clash_distance, skip_last):
    res_i = list(get_chain_residue_indices(pose, chain_i, residues_per_chain))
    res_j = list(get_chain_residue_indices(pose, chain_j, residues_per_chain))

    # Skip last N residues
    if len(res_i) > skip_last:
        res_i = res_i[:-skip_last]
    if len(res_j) > skip_last:
        res_j = res_j[:-skip_last]
    
    backbone_atoms = {"N", "CA", "C", "O"}

    coords_i, atom_indices_i = [], []
    for r1 in res_i:
        res1 = pose.residue(r1)
        for a1 in range(1, res1.natoms() + 1):
            atom_name = res1.atom_name(a1).strip()
            if atom_name in backbone_atoms:
                vec = res1.xyz(a1)
                if is_valid_xyz(vec):
                    coords_i.append(vec)
                    atom_indices_i.append((r1, a1))

    coords_j, atom_indices_j = [], []
    for r2 in res_j:
        res2 = pose.residue(r2)
        for a2 in range(1, res2.natoms() + 1):
            atom_name = res2.atom_name(a2).strip()
            if atom_name in backbone_atoms:
                vec = res2.xyz(a2)
                if is_valid_xyz(vec):
                    coords_j.append(vec)
                    atom_indices_j.append((r2, a2))

    coords_i_np = np.array([[v.x, v.y, v.z] for v in coords_i])
    coords_j_np = np.array([[v.x, v.y, v.z] for v in coords_j])

    logger.info(f"coords_i_np: {len(coords_i_np)}, atom_indices_i: {len(atom_indices_i)}")
    logger.info(f"coords_j_np: {len(coords_j_np)}, atom_indices_j: {len(atom_indices_j)}")

    if coords_i_np.shape[0] == 0 or coords_j_np.shape[0] == 0:
        raise ValueError("One of the coordinate arrays is empty. Cannot build KDTree.")

    tree_j = cKDTree(coords_j_np)

    # output_dir = Path("/home/d12-studenti/Kiyan/prosculpt/test_plugins/closest_pairs")
    # output_dir.mkdir(parents=True, exist_ok=True)

    # pose_name = pose.pdb_info().name() or "unknown"
    # pose_name_clean = Path(pose_name).stem.replace(" ", "_").replace("/", "_")
    # r1, a1 = atom_indices_i[idx_i]
    # r2, a2 = atom_indices_j[idx_j]
    # clash_tag = output_dir / f"{pose_name_clean}_res{r1}_{r2}.pdb"

    # pose.dump_pdb(str(clash_tag))
    # logger.info(f"Clash found. Dumped structure to {clash_tag}")

    for idx_i, coord in enumerate(coords_i_np):
        dist, idx_j = tree_j.query(coord, k=1)

        if dist < clash_distance:
            # Logging clash detail
            r1, a1 = atom_indices_i[idx_i]
            r2, a2 = atom_indices_j[idx_j]
            name1 = pose.residue(r1).atom_name(a1).strip()
            name2 = pose.residue(r2).atom_name(a2).strip()
            logger.warning(f"Clash detected: {name1} (res {r1}) ↔ {name2} (res {r2}), distance: {dist:.2f} Å")

            return True, dist  # Early return on first clash
    
    # # If no clash detected, log the closest pair
    # if (
    #     min_idx_i is not None and min_idx_j is not None and
    #     min_idx_i < len(atom_indices_i) and min_idx_j < len(atom_indices_j)
    # ):
    #     r1, a1 = atom_indices_i[min_idx_i]
    #     r2, a2 = atom_indices_j[min_idx_j]
    #     name1 = pose.residue(r1).atom_name(a1).strip()
    #     name2 = pose.residue(r2).atom_name(a2).strip()
    #     logger.info(f"No clash detected. Closest atoms: {name1} (res {r1}) ↔ {name2} (res {r2}), distance: {min_dist_overall:.2f} Å")

    #      # Save a visualization PDB for inspection
    #     output_dir = Path("/home/d12-studenti/Kiyan/prosculpt/test_plugins/closest_pairs")
    #     output_dir.mkdir(parents=True, exist_ok=True)

    #     pose_name = pose.pdb_info().name() or "unknown"
    #     pose_name_clean = Path(pose_name).stem.replace(" ", "_").replace("/", "_")
    #     close_tag = output_dir / f"closest_{pose_name_clean}_res{r1}_{r2}_dist{min_dist_overall:.2f}.pdb"

    #     pose.dump_pdb(str(close_tag))
    #     logger.info(f"Dumped closest-pair PDB to {close_tag}")
    
    # else:
    #     logger.warning("Could not log closest atoms due to index mismatch.")

    return False

def has_any_inter_chain_clashes(pose, residues_per_chain, clash_distance, skip_last, num_chains=8):
    for i in range(num_chains):
        j = (i + 1) % num_chains
        clash = has_chain_pair_clashes(
            pose, i, j, residues_per_chain, clash_distance, skip_last=skip_last
        )
        if clash:
            return True
    return False, float('inf')

# Plugin entry point
def filter_backbone(pdb: str | Path, kwargs: dict):
    ensure_pyrosetta_initialized()

    base_path = kwargs["base_pdb"]
    residues_to_align = kwargs.get("residues_to_align", 85)
    clash_distance = kwargs.get("clash_distance", 0.3)
    num_chains = kwargs.get("num_chains", 8)
    skip_last = kwargs.get("skip_last_residues", 0)

    logger.info(f"Running clash_check on: {pdb} with settings: base_pdb={base_path}, residues_to_align={residues_to_align}, clash_distance={clash_distance}, num_chains={num_chains}, skip_last_residues={skip_last}")

    design = pose_from_file(str(pdb))
    residues_per_chain = len(get_chain_residues(design, 'A'))

    combined = assemble_c8_design(base_path, str(pdb), residues_to_align=residues_to_align)

    has_clash = has_any_inter_chain_clashes(
        combined,
        residues_per_chain=residues_per_chain,
        clash_distance=clash_distance,
        skip_last=skip_last,
        num_chains=num_chains
    )

    return not has_clash
