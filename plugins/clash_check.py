from pyrosetta import *
from pyrosetta import Pose, pose_from_file
from pyrosetta.rosetta.utility import vector1_bool
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.protocols.toolbox.pose_manipulation import superimpose_pose_on_subset_CA
from pyrosetta.rosetta.utility import vector1_unsigned_long
from pathlib import Path
import logging

# from plugin_utils import ensure_pyrosetta_initialized

logger = logging.getLogger(__name__)

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

def assemble_c8_design(base_pdb: str, design_pdb: str, residues_to_align: int):
    base = pose_from_file(base_pdb)
    design = pose_from_file(design_pdb)

    design_chain_A_residues = get_chain_residues(design, 'A')
    if len(design_chain_A_residues) < residues_to_align:
        raise ValueError("Design chain A is too short for alignment.")

    design_align_residues = design_chain_A_residues[-residues_to_align:]
    num_res_base_subunit = len(get_chain_residues(base, 'A'))
    combined = None

    for i in range(2):
        base_start = 1 + i * num_res_base_subunit
        base_align_vec = make_vector1_range(
            base_start + num_res_base_subunit - residues_to_align,
            base_start + num_res_base_subunit - 1
        )
        aligned_copy = copy_pose(design)
        offset = -1 * ((base_start + num_res_base_subunit - residues_to_align) - design_align_residues[0])
        superimpose_pose_on_subset_CA(aligned_copy, base, base_align_vec, offset=offset)
        aligned_copy.pdb_info().set_chains(chr(ord('A') + i))
        if combined is None:
            combined = copy_pose(aligned_copy)
        else:
            combined.append_pose_by_jump(aligned_copy, 1)

    return combined

def get_clash_xml():
    return """
<ROSETTASCRIPTS>
    
        <SCOREFXNS>
            <ScoreFunction name="sfxn" weights="ref2015" />
        </SCOREFXNS>

        <RESIDUE_SELECTORS>
            <ScoreTermValueBased name="clashing_res" score_type="fa_rep" score_fxn="sfxn" lower_threshold="10" upper_threshold="99999" />
            <Chain name="only_chain_A" chains="A" />
            <Chain name="only_chain_B" chains="B" />
            <And name="clash_in_A" selectors="clashing_res,only_chain_A" />
            <And name="clash_in_B" selectors="clashing_res,only_chain_B" />
            <Or name="clash_in_A_and_B" selectors="clash_in_A, clash_in_B" />
        </RESIDUE_SELECTORS>

        <FILTERS>
            <ResidueCount name="clash_count_A" residue_selector="clash_in_A" />
            <ResidueCount name="clash_count_B" residue_selector="clash_in_B" />
            <ResidueCount name="clash_count_A_and_B" residue_selector="clash_in_A_and_B" />
        </FILTERS>

        <PROTOCOLS>
            <Add filter_name="clash_count_A" />
            <Add filter_name="clash_count_B" />
            <Add filter_name="clash_count_A_and_B" />
        </PROTOCOLS>

    </ROSETTASCRIPTS>
"""

# Plugin entry point
def filter_backbone(pdb: str | Path, kwargs: dict):
    """
    Checks for inter-chain clashes between the designed structure and the reference base structure.

    Parameters:
        pdb (str | Path): Path to the PDB file to be evaluated.

        kwargs (dict): Plugin-specific settings from the YAML configuration, expected to include:
            - base_pdb (str): Path to the base/reference PDB file used for alignment and clash comparison.
            - residues_to_align (int): Number of C-terminal residues in chain A to use for alignment.
            - skip_last_residues (int): Number of residues at the C-terminus to exclude from clash checking.
    """

    # ensure_pyrosetta_initialized()

    base_path = kwargs["base_pdb"]
    residues_to_align = kwargs.get("residues_to_align", 85)
    skip_last = kwargs.get("skip_last_residues", 200)

    logger.info(f"Running RosettaScripts clash_check on: {pdb} (base: {base_path}, align: {residues_to_align}, skip_last: {skip_last})")

    design = pose_from_file(str(pdb))
    combined = assemble_c8_design(base_path, str(pdb), residues_to_align=residues_to_align)

    xml_obj = XmlObjects.create_from_string(get_clash_xml())
    clash_A_sel = xml_obj.get_residue_selector("clash_in_A")
    clash_both_sel = xml_obj.get_residue_selector("clash_in_A_and_B")

    subset_A = clash_both_sel.apply(design)
    subset_both = clash_A_sel.apply(combined)

    limit = design.size() - skip_last
    subset_A_filtered = vector1_bool(design.size())
    for i in range(1, limit + 1):
        if subset_A[i]:
            subset_A_filtered[i] = True

    design_len = design.size()
    limit1 = design_len - skip_last

    subset_both_filtered = vector1_bool(combined.size())
    for i in range(1, combined.size() + 1):
        if subset_both[i]:
            # Keep residue if it's not in last 200 of either subunit
            if (i <= limit1):
                subset_both_filtered[i] = True


    logger.info("Residues with clashes in chain A:")
    logger.info([i for i, is_selected in enumerate(subset_A, start=1) if is_selected])

    logger.info("Residues with clashes in both chains:")
    logger.info([i for i, is_selected in enumerate(subset_both, start=1) if is_selected])

    logger.info(f"Filtered residues with clashes in chain A (before last {skip_last}):")
    logger.info([i for i, selected in enumerate(subset_A_filtered, start=1) if selected])

    logger.info(f"Filtered residues with clashes in both chains (excluding last {skip_last} of each):")
    logger.info([i for i, selected in enumerate(subset_both_filtered, start=1) if selected])

    for i in range(1, design.size() + 1):
        residue = design.residue(i)

        if not residue.is_protein() or not residue.has("CA"):
            continue

        try:
            ca_idx = residue.atom_index("CA")

            if subset_both_filtered[i] and not subset_A_filtered[i]:
                combined.pdb_info().bfactor(i, ca_idx, 15.0)  # Orange/Yellow or custom color

        except Exception as e:
            logger.warning(f"Error tagging residue {i} ({residue.name()}): {e}")
    

    # Dump the combined pose to a PDB file
    # output_dir = Path("/home/d12-studenti/Kiyan/prosculpt/test_plugins/assembled_rollers")
    # output_dir.mkdir(parents=True, exist_ok=True)
    # output_path = output_dir / f"assembled_{Path(pdb).stem}.pdb"
    # combined.dump_pdb(str(output_path))
    # logger.info(f"Dumped assembled pose to: {output_path}")


    count_A = sum(subset_A_filtered)
    count_both = sum(subset_both_filtered)
    inter_chain_clashes = count_both - count_A

    logger.info(f"Detected {inter_chain_clashes} estimated inter-chain clashes")

    return inter_chain_clashes == 0
