from pathlib import Path
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.select.residue_selector import (
    ChainSelector, AndResidueSelector, NotResidueSelector,
    NeighborhoodResidueSelector, ResidueIndexSelector
)
from pyrosetta.rosetta.protocols.moves import DsspMover

# from plugin_utils import ensure_pyrosetta_initialized

import logging

logger = logging.getLogger("prosculpt.plugins")

def get_first_helix_range(pose, chain_id="A"):
    DsspMover().apply(pose)
    dssp = pose.secstruct()

    helix_start = None
    helix_end = None

    for i in range(1, pose.size() + 1):
        if pose.pdb_info().chain(i) != chain_id:
            continue
        if dssp[i - 1] == 'H':
            if helix_start is None:
                helix_start = i
            helix_end = i
        elif helix_start is not None:
            break

    if helix_start is None:
        raise ValueError(f"No helix found in chain {chain_id}")

    return f"{helix_start}-{helix_end}"


def filter_backbone(pdb: str | Path, kwargs: dict):
    """
    Filter that returns True if the number of neighbor residues is within the specified range.
    """

    # ensure_pyrosetta_initialized()
    logger.info(f"Running num_neighbors on: {pdb}")

    pose = pose_from_pdb(str(pdb))
    pose.update_residue_neighbors()

    # Get helix range
    helix_range = get_first_helix_range(pose, chain_id="A")

    # === Define residue selectors ===
    seed_index = ResidueIndexSelector(helix_range)
    seed_chain = ChainSelector("A")
    seed_selector = AndResidueSelector(seed_index, seed_chain)

    scaffold_selector = ChainSelector("B")
    non_scaffold_selector = NotResidueSelector(scaffold_selector)

    neighbor_selector = NeighborhoodResidueSelector()
    neighbor_selector.set_focus_selector(seed_selector)

    neighbor_distance = kwargs.get("neighbor_distance", 8)
    if neighbor_distance is None:
        raise ValueError("Missing required parameter: 'neighbor_distance'")
    neighbor_selector.set_distance(neighbor_distance)
    
    neighbor_selector.set_include_focus_in_subset(False)

    final_neighbor_selector = AndResidueSelector(neighbor_selector, non_scaffold_selector)
    neighbor_subset = final_neighbor_selector.apply(pose)

    num_neighbors = sum(neighbor_subset)

    # Get thresholds
    min_neighbors = kwargs.get("min_num_neighbors", 0)
    max_neighbors = kwargs.get("max_num_neighbors", float('inf'))
    

    return min_neighbors <= num_neighbors <= max_neighbors
