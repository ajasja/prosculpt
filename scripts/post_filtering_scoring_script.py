"""Example post_filtering_scoring_script (see the `filtering:` config section
and run_post_filtering_scoring() in prosculpt_run.py). Runs Rosetta's
InterfaceAnalyzerMover and ShapeComplementarity filter for one or more
interfaces against prosculpt's already-filtered subset of models, in
parallel across the job's available CPUs. An interface can be chain-to-chain
("A_B") or group-to-group ("AB_CD").

Does not recompute Rg/charge/SAP - those are already computed by the
mandatory scoring_script.py pass earlier in the pipeline.

Always invoked the same way (`your_script.py <filtered_output.csv>
--output_dir <dir>`, plus any extra flags from a job's own
post_filtering_scoring_script_arguments). USE_BACKBONE_MINIMIZATION and
CHAIN_PAIRS below are the file-level defaults, overridable per job via
--chain_pairs/--use_backbone_minimization.
"""

from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pyrosetta import init, pose_from_pdb, create_score_function
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.kinematics import MoveMap
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
import pandas as pd
import argparse
import itertools
import os
import time

# ============================================================================
# User-configurable settings. Edit these two values directly.
# ============================================================================

# Whether to run a backbone minimization step before interface analysis
# (side chain -> backbone -> side chain, Cartesian) or side-chain-only.
USE_BACKBONE_MINIMIZATION = True

# Which interfaces to score ddG / shape-complementarity for. An interface is
# written as the two sides separated by an underscore:
#   - "all"        every pairwise combination of chains present in each model
#   - "A_B"        chain A against chain B
#   - "A_B,A_C"    a comma-separated list of interfaces
#   - "AB_CD"      chains A+B together against chains C+D together
#   - "A_BCD"      chain A against the other three as one group
# An interface is silently skipped for a model if any of its chains isn't
# present in it. "all" only enumerates pairs - name group interfaces
# explicitly, since there are too many possible groupings to enumerate.
CHAIN_PAIRS = "A_B"

def _init_worker():
    """Runs once per worker process to initialize PyRosetta - its C++ state
    doesn't survive being forked/pickled, so each worker needs its own init()."""
    init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy -corrections::beta_nov16 true")

def get_protein_chain_ids(structure):
    """
    Returns the list of chain IDs (in file order) that contain at least one
    amino acid residue, i.e. skips purely hetero/water chains.

    Parameters:
        structure (Bio.PDB.Structure.Structure): A parsed structure (first model is used).

    Returns:
        list[str]: Chain IDs that contain protein residues.
    """
    model = structure[0]
    return [chain.id for chain in model if any(is_aa(res, standard=False) for res in chain)]

def parse_chain_pairs(pairs_arg):
    """
    Parses the CHAIN_PAIRS setting into either "all" or a list of
    (group1, group2) interfaces, each side a list of chain letters:
    "A_B"        -> [(["A"], ["B"])]
    "A_B,A_C"    -> [(["A"], ["B"]), (["A"], ["C"])]
    "AB_CD"      -> [(["A", "B"], ["C", "D"])]

    Parameters:
        pairs_arg (str): CHAIN_PAIRS, e.g. "all", "A_B", "A_B,A_C" or "AB_CD".

    Returns:
        "all" or list[tuple[list[str], list[str]]]: The parsed interfaces, in the order given.
    """
    if pairs_arg.strip().lower() == "all":
        return "all"

    interfaces = []
    for token in pairs_arg.split(","):
        token = token.strip()
        if not token:
            continue
        if "_" not in token:
            hint = f" Did you mean {'_'.join(token)!r}?" if len(token) == 2 else ""
            raise ValueError(
                f"Invalid interface {token!r}: the two sides of an interface are separated by an "
                f"underscore, e.g. 'A_B' or 'AB_CD'.{hint}"
            )
        left, _, right = token.partition("_")
        group1, group2 = list(left.strip()), list(right.strip())
        if not group1 or not group2:
            raise ValueError(f"Invalid interface {token!r}: both sides need at least one chain.")
        shared = set(group1) & set(group2)
        if shared:
            raise ValueError(
                f"Invalid interface {token!r}: chain(s) {sorted(shared)} appear on both sides."
            )
        for side in (group1, group2):
            if len(set(side)) != len(side):
                raise ValueError(f"Invalid interface {token!r}: a chain is repeated on one side.")
        interfaces.append((group1, group2))
    return interfaces

def minimize_pose(pose, use_backbone_minimization=True):
    """
    Minimizes a pose in place: side chain -> backbone -> side chain
    (Cartesian) if backbone minimization is requested, otherwise
    side-chain-only (non-Cartesian).

    Parameters:
        pose (pyrosetta.Pose): The pose to minimize, modified in place.
        use_backbone_minimization (bool): Whether to also minimize the backbone.

    Returns:
        str: The name of the score function weights used for minimization
             (and that should be used for the subsequent interface analysis).
    """
    scorefxn_name = "beta_nov16_cart" if use_backbone_minimization else "beta_nov16"
    scorefxn = create_score_function(scorefxn_name)

    sc_movemap = MoveMap()
    sc_movemap.set_chi(True)
    sc_movemap.set_bb(False)
    minimize_sc = MinMover(sc_movemap, scorefxn, "lbfgs_armijo_nonmonotone", 0.01, True)
    minimize_sc.cartesian(use_backbone_minimization)
    minimize_sc.apply(pose)

    if use_backbone_minimization:
        bb_movemap = MoveMap()
        bb_movemap.set_chi(True)
        bb_movemap.set_bb(True)
        minimize_bb = MinMover(bb_movemap, scorefxn, "lbfgs_armijo_nonmonotone", 0.01, True)
        minimize_bb.cartesian(True)
        minimize_bb.apply(pose)
        minimize_sc.apply(pose)

    return scorefxn_name

def interface_label(group1, group2):
    """
    The column-name prefix for one interface, e.g. "A_B" or "AB_CD".

    Parameters:
        group1 (list[str]): Chain letters on one side.
        group2 (list[str]): Chain letters on the other side.

    Returns:
        str: The two sides joined by an underscore, as configured.
    """
    return f"{''.join(group1)}_{''.join(group2)}"

def chain_letters_in_order(pose):
    """
    The pose's PDB chain letters, deduplicated, in the order residues appear.

    SwitchChainOrder renames the chains it keeps (the first kept chain becomes
    A, the second B, ...), so the letters an extracted pose actually has are
    not the ones that were asked for. Reading them back off the pose keeps the
    residue selectors below correct either way.

    Parameters:
        pose (pyrosetta.Pose): Any pose with PDBInfo.

    Returns:
        list[str]: Chain letters, first-seen order.
    """
    info = pose.pdb_info()
    letters = []
    for resi in range(1, pose.total_residue() + 1):
        letter = info.chain(resi)
        if letter not in letters:
            letters.append(letter)
    return letters

def calculate_interface_metrics(pose, group1, group2, scorefxn_name, use_backbone_minimization):
    """
    Extracts the two groups of chains from an already-minimized pose into
    their own pose and runs InterfaceAnalyzerMover plus the
    ShapeComplementarity filter across the group1/group2 boundary.

    InterfaceAnalyzerMover rather than the ddG mover because the latter
    separates the pose across a single jump, which cannot express a
    group-vs-group boundary - it would silently report one chain against all
    the others. IA takes the two sides explicitly.

    Parameters:
        pose (pyrosetta.Pose): The full, already-minimized pose (not modified).
        group1 (list[str]): Chain letters on one side of the interface.
        group2 (list[str]): Chain letters on the other side.
        scorefxn_name (str): Score function weights to use (matches minimize_pose's).
        use_backbone_minimization (bool): Whether the interface minimization
            should be Cartesian (must match scorefxn_name).

    Returns:
        dict: Scores prefixed with the interface label, e.g.
              {"AB_CD_dG_separated": ..., "AB_CD_sc_value": ..., "AB_CD_sc2": ...}
              for "AB_CD", or {"A_B_dG_separated": ..., ...} for "A_B".
    """
    prefix = interface_label(group1, group2)
    cartesian = "true" if use_backbone_minimization else "false"

    interface_pose = pose.clone()
    kept_chains = group1 + group2
    extract = XmlObjects.static_get_mover(
        f'<SwitchChainOrder name="take_interface" chain_name="{",".join(kept_chains)}"/>'
    )
    extract.apply(interface_pose)

    # Post-extraction letters, split back into the two sides by group size -
    # the order chains were asked for above is the order they end up in.
    letters = chain_letters_in_order(interface_pose)
    if len(letters) != len(kept_chains):
        raise RuntimeError(
            f"Expected {len(kept_chains)} chain(s) after extracting {prefix}, got {letters}."
        )
    side1_chains = letters[: len(group1)]
    side2_chains = letters[len(group1) :]
    side1 = ",".join(side1_chains)
    side2 = ",".join(side2_chains)
    # InterfaceAnalyzerMover's dock_chains-style interface string, e.g. "AB_CD".
    ia_interface = f"{''.join(side1_chains)}_{''.join(side2_chains)}"

    xml = f"""
    <ROSETTASCRIPTS>
        <SCOREFXNS>
            <ScoreFunction name="sfxn_clean" weights="{scorefxn_name}" />
        </SCOREFXNS>
        <RESIDUE_SELECTORS>
            <Chain name="side1" chains="{side1}" />
            <Chain name="side2" chains="{side2}" />
        </RESIDUE_SELECTORS>
        <FILTERS>
            <ShapeComplementarity name="sc2" min_sc="0.0" verbose="0" quick="0"
                residue_selector1="side1" residue_selector2="side2"
                write_int_area="1" write_median_dist="1" confidence="0" />
        </FILTERS>
        <MOVERS>
            <MinMover name="minimize_sc_all" scorefxn="sfxn_clean" bb="0" chi="1" cartesian="{cartesian}" />
            <InterfaceAnalyzerMover name="analyze_interface" scorefxn="sfxn_clean"
                interface="{ia_interface}"
                packstat="0" interface_sc="1" pack_input="0" pack_separated="1" />
        </MOVERS>
        <PROTOCOLS>
            <Add mover="minimize_sc_all" />
            <Add mover="analyze_interface" />
            <Add filter="sc2" />
        </PROTOCOLS>
    </ROSETTASCRIPTS>
    """
    protocol = XmlObjects.create_from_string(xml).get_mover("ParsedProtocol")
    protocol.apply(interface_pose)

    return {f"{prefix}_{key}": value for key, value in dict(interface_pose.scores).items()}
def score_one_model(pdb, chain_pairs, use_backbone_minimization=USE_BACKBONE_MINIMIZATION):
    """
    Computes, for each configured interface present in the model,
    ddG/shape-complementarity for a single PDB model.

    Runs inside a worker process (see _init_worker) - its arguments and
    return value must be picklable. chain_pairs and use_backbone_minimization
    are passed explicitly (via main()'s partial()) rather than read from the
    module-level globals, since a 'spawn'ed worker re-imports the module and
    wouldn't see a --chain_pairs/--use_backbone_minimization override.

    Parameters:
        pdb (str): Path to the model's PDB file.
        chain_pairs ("all" or list[tuple[list[str], list[str]]]): Interfaces to
            score, as returned by parse_chain_pairs. "all" scores every
            pairwise combination of chains present in this model (group
            interfaces have to be named explicitly - there are too many
            possible groupings to enumerate usefully).
        use_backbone_minimization (bool): Whether to also minimize the
            backbone before interface analysis (see minimize_pose()).

    Returns:
        dict: One row's worth of scores, ready for pd.DataFrame(..., index=[0]).
    """
    start_time = time.time()

    bio_parser = PDBParser(QUIET=True)
    structure = bio_parser.get_structure('protein', pdb)

    pose = pose_from_pdb(pdb)

    dictionary = {'model_path': pdb}

    protein_chain_ids = get_protein_chain_ids(structure)
    if chain_pairs == "all":
        interfaces_to_run = [([c1], [c2]) for c1, c2 in itertools.combinations(protein_chain_ids, 2)]
    else:
        interfaces_to_run = []
        for group1, group2 in chain_pairs:
            missing = [c for c in group1 + group2 if c not in protein_chain_ids]
            if missing:
                print(
                    f"Skipping interface {interface_label(group1, group2)} for {pdb}: "
                    f"chain(s) {missing} not found (chains present: {protein_chain_ids})."
                )
                continue
            interfaces_to_run.append((group1, group2))

    if interfaces_to_run:
        scorefxn_name = minimize_pose(pose, use_backbone_minimization=use_backbone_minimization)
        for group1, group2 in interfaces_to_run:
            dictionary.update(
                calculate_interface_metrics(pose, group1, group2, scorefxn_name, use_backbone_minimization)
            )

    # Timed inside the worker process, since only it knows this model's own wall-clock cost.
    elapsed = time.time() - start_time
    print(f"Scored {os.path.basename(pdb)} in {elapsed:.1f}s ({len(interfaces_to_run)} interface(s)).", flush=True)

    return dictionary

def _parse_bool_arg(value):
    """
    Parses a command-line boolean ("true"/"false", "1"/"0", "yes"/"no", any
    case) - argparse's own type=bool is a trap, since bool("false") is True.
    """
    normalized = value.strip().lower()
    if normalized in ("true", "1", "yes"):
        return True
    if normalized in ("false", "0", "no"):
        return False
    raise argparse.ArgumentTypeError(f"Expected a boolean (true/false), got {value!r}.")


def main():
    parser = argparse.ArgumentParser(description='Calculate (per CHAIN_PAIRS/USE_BACKBONE_MINIMIZATION above, or --chain_pairs/--use_backbone_minimization) per-interface ddG/shape-complementarity for a list of PDB models.')
    parser.add_argument('input_csv', type=str, help='Path to input CSV file')
    parser.add_argument('--output_dir', type=str, help='Path to output directory. If not provided, output files will be written to the same directory as the input CSV.')
    parser.add_argument('--num_workers', type=int, default=None,
        help='Number of models to score in parallel. Defaults to the SLURM_CPUS_PER_TASK allocation if set, otherwise the number of locally visible CPUs.')
    parser.add_argument('--chain_pairs', type=str, default=None,
        help='Overrides the CHAIN_PAIRS setting above for this run only. An interface is the two sides separated by an underscore: "A_B" (chain against chain), "AB_CD" (chains A+B against C+D), a comma-separated list like "A_B,AB_CD", or "all" for every pairwise combination present. Defaults to CHAIN_PAIRS if not given.')
    parser.add_argument('--use_backbone_minimization', type=_parse_bool_arg, default=None,
        help='Overrides the USE_BACKBONE_MINIMIZATION setting above for this run only (true/false). Defaults to USE_BACKBONE_MINIMIZATION if not given.')
    args = parser.parse_args()
    chain_pairs = parse_chain_pairs(args.chain_pairs if args.chain_pairs is not None else CHAIN_PAIRS)
    use_backbone_minimization = (
        args.use_backbone_minimization if args.use_backbone_minimization is not None else USE_BACKBONE_MINIMIZATION
    )

    data = pd.read_csv(args.input_csv)

    if args.output_dir:
        output = args.output_dir
    else:
        output = os.path.dirname(args.input_csv)

    num_workers = args.num_workers or int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
    num_workers = max(1, num_workers)
    if len(data) > 0:
        num_workers = min(num_workers, len(data))

    path_csv = os.path.join(output, "post_filtering_scores.csv")

    start_time = time.time()
    worker_fn = partial(score_one_model, chain_pairs=chain_pairs, use_backbone_minimization=use_backbone_minimization)
    with ProcessPoolExecutor(max_workers=num_workers, initializer=_init_worker) as executor:
        for dictionary in executor.map(worker_fn, data['model_path']):
            df = pd.DataFrame(dictionary, index=[0])
            df.to_csv(path_csv, mode='a', header=not os.path.exists(path_csv), index=False)
    elapsed = time.time() - start_time
    per_model = elapsed / len(data) if len(data) else 0.0
    print(
        f"Post-filtering scoring: scored {len(data)} model(s) using {num_workers} worker(s) "
        f"in {elapsed:.1f}s total ({elapsed / 60:.1f} min) - {per_model:.1f}s/model on average.",
        flush=True,
    )

if __name__ == "__main__":
    main()
