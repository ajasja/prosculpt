from omegaconf import DictConfig, OmegaConf
import hydra
from hydra.core.hydra_config import HydraConfig
import os
from hydra.utils import get_original_cwd, to_absolute_path
import logging
import prosculpt
import glob
import re
import shutil
import pathlib
from pathlib import Path
from omegaconf import open_dict
from Bio import SeqIO
from pyrosetta import *

# from plugins.plugin import filter_backbones_after_rfdiff
from Bio.PDB import PDBParser

import importlib.util
import logging
import shutil
from pathlib import Path
from typing import List
import pprint


log = logging.getLogger(__name__)

init("-mute all")


def run_and_log(command, log_func=log.info, dry_run=False, cfg=None):
    """Runs a command using os.system and also logs the command before running using log.info"""
    if log_func:
        log_func(command)
    if not dry_run:
        stat = os.system(command)
        wife = os.WIFEXITED(stat)
        exitCode = os.waitstatus_to_exitcode(stat)
        log.info(
            f"Command exited with status {stat} and WIFEXITED {wife}. Exit code: {exitCode}"
        )
        if exitCode != 0:
            log.error(
                "There was an error running the command. We consider it fatal to prevent any file loss. Check the logs and contact the developer."
            )
            dodatek = ""

            auto_restart_count = cfg.get("auto_restart_count", 0)
            auto_restart = cfg.get("auto_restart", 0)
            log.info(
                f"auto_restart_count: {auto_restart_count}, auto_restart: {auto_restart}"
            )
            if auto_restart_count < auto_restart:  # if cfg has autoRestart > 0
                log.error("As per config, Prosculpt will be restarted automatically.")
                with open_dict(cfg):
                    cfg.auto_restart_count = auto_restart_count + 1
                dodatek = f"However, Prosculpt was autoRestarted {cfg.get('auto_restart_count', 0)} out of {cfg.get('auto_restart', 0)} times (as per config) after encountering the crash."
                prosculptApp(cfg)
            raise Exception(f"Command exited with exit code {exitCode}\n\n{dodatek}")
        # TODO: make `cfg` global and manipulate it everywhere. To avoid the need of passing `cfg=cfg` to every run_and_log call.


scripts_folder = pathlib.Path(__file__).resolve().parent / "scripts"


def general_config_prep(cfg):
    log.info("Running generalPrep")
    os.makedirs(cfg.output_dir, exist_ok=True)

    # We want to have all the params in the cfg struct. Thus, use we open_dict to be able to write new data.
    with open_dict(cfg):

        cfg.output_dir = str(
            pathlib.Path(cfg.output_dir).resolve()
        )  # Output dir MUST be absolute, otherwise ProtMPPN complains about missing seq_chain (update: it probably doesn't matter)
        cfg.rfdiff_out_dir = os.path.join(cfg.output_dir, "1_rfdiff")
        cfg.mpnn_out_dir = os.path.join(cfg.output_dir, "2_mpnn")
        cfg.af2_out_dir = os.path.join(cfg.output_dir, "3_models")
        cfg.rfdiff_out_path = os.path.join(
            cfg.rfdiff_out_dir, ""
        )  # Appending empty string "" results in rfdiff_out_path ending with a correct directory separator (/ or \)

        cfg.path_for_parsed_chains = os.path.join(cfg.mpnn_out_dir, "parsed_pdbs.jsonl")
        cfg.path_for_assigned_chains = os.path.join(
            cfg.mpnn_out_dir, "assigned_pdbs.jsonl"
        )
        cfg.path_for_fixed_positions = os.path.join(
            cfg.mpnn_out_dir, "fixed_pdbs.jsonl"
        )
        cfg.path_for_tied_positions = os.path.join(cfg.mpnn_out_dir, "tied_pdbs.jsonl")

        cfg.fasta_dir = os.path.join(cfg.mpnn_out_dir, "seqs")
        cfg.rfdiff_pdb = os.path.join(cfg.rfdiff_out_path, "_0.pdb")

        cfg.chains_to_design = cfg.get(
            "chains_to_design", None
        )  # If not specified, will be determined from contig later
        cfg.use_a3m = cfg.get("use_a3m", False)  # whether to use a3m files
        cfg.prediction_model = cfg.get(
            "prediction_model", "Colabfold"
        )  # Options are Colabfold or Boltz2
        print(f"Prediction model: {cfg.prediction_model}")

        if cfg.get("skipRfDiff", False):
            # We only need to redesign the chains specified in designable_residues
            cfg.chains_to_design = " ".join(
                sorted({_[0] for _ in cfg.designable_residues})
            )
            log.info(
                f"Skipping RFdiff, only redesigning chains specified in designable_residues: {cfg.chains_to_design}"
            )

        if "inference" not in cfg:
            cfg.inference = {}

        if "symmetry" not in cfg.inference:
            cfg.inference.symmetry = None

        if "omit_AAs" not in cfg:
            cfg.omit_AAs = "X"  # This is the default also in proteinMPNN

        # I suggest the following: count("/0", contig) -> chains_to_design = " ".join(chain_letters[:count]), unless specified (in run.yaml, it should be null, None or sth similar)
        chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"  # What happens after 26 chains? RfDiff only supports up to 26 chains: https://github.com/RosettaCommons/RFdiffusion/blob/ba8446eae0fb80c121829a67d3464772cc827f01/rfdiffusion/contigs.py#L40C29-L40C55
        if cfg.chains_to_design == None:  # TODO this will likely break for sym mode

            if cfg.inference.symmetry != None:  # if we are doing symmetry:
                breaks = int(cfg.inference.symmetry[1:])
            else:
                breaks = cfg.contig.count("/0 ") + 1
            cfg.chains_to_design = " ".join(chain_letters[:breaks])
            log.info(
                f"Chains to design (according to contig chain breaks): {cfg.chains_to_design}"
            )

    for directory in [cfg.rfdiff_out_dir, cfg.mpnn_out_dir, cfg.af2_out_dir]:
        os.makedirs(directory, exist_ok=True)
        log.info(f"Made directory {directory}")

    # No need to return; cfg is mutable.


def run_rfdiff(cfg):
    """
    RFdiffusion will generate new protein structures according to the contig specified
    INPUT: starting pdb, contig, number of structures to design
    OUTPUT: generated structure (pdb file), metadata associated with specific run for each generated structure (trb format)

    TRB file contains useful information.
    In this script data from con_hal_pdb_idx/complex_con_hal_pdb_idx and 'complex_con_ref_idx0' are used in helper functions
    See RFdiffusion git for details.
    """
    log.info("***************Running runRFdiff***************")

    # Check if we already have cfg.num_designs_rfdiff .pdb and *.trb files in cfg.rfdiff_out_path
    if (
        len(glob.glob(os.path.join(cfg.rfdiff_out_path, "*.pdb")))
        == cfg.num_designs_rfdiff
    ):
        log.info(
            f"Found {cfg.num_designs_rfdiff} .pdb and .trb files in {cfg.rfdiff_out_path}. Skipping RFdiff."
        )
        if (
            len(glob.glob(os.path.join(cfg.rfdiff_out_path, "*.trb")))
            != cfg.num_designs_rfdiff
        ):
            log.critical(
                f"Found {len(glob.glob(os.path.join(cfg.rfdiff_out_path, '*.trb')))} trb files in {cfg.rfdiff_out_path}, expected {cfg.num_designs_rfdiff}. \n Most likely, your previous simulation crashed while RfDiff was writing output files. Please, manually remove the .pdb file which does not have an associated .trb file in the {cfg.rfdiff_out_path} directory, then run the computation again."
            )
            raise Exception(
                f"Number of RfDiff .pdb files does not match RfDiff .trb files!"
            )
        # endif
        log.info("***************Skipping RFdiffusion altogether***************")
        return

    if r"{OUTPUT_PATH}" in cfg.rfdiff_run_command:
        rfdiff_run_command = (cfg.rfdiff_run_command).replace(
            r"{OUTPUT_PATH}", cfg.output_dir
        )

        out_path = "/output"
        rfdiff_out_path = "/output/1_rfdiff/"
    else:
        rfdiff_run_command = cfg.rfdiff_run_command
        rfdiff_out_path = cfg.rfdiff_out_path
        out_path = cfg.output_dir

    # os.makedirs(rfdiff_out_path, exist_ok=True)
    os.makedirs(cfg.output_dir + "/schedules", exist_ok=True)

    rfdiff_cmd_str = f"{rfdiff_run_command} \
          inference.output_prefix={rfdiff_out_path} \
          'contigmap.contigs={cfg.contig}' \
          inference.num_designs={cfg.num_designs_rfdiff} \
          -cn prosculpt2rfdiff.yaml -cd {out_path}"
    run_and_log(rfdiff_cmd_str, cfg=cfg)

    log.info("***************After running RFdiffusion***************")


def rechain_rfdiff_pdbs(cfg):
    """
    RFdiffusion joins chain sequences together
    For contig "[A1-30/4-6/C1-30/0 D1-30/0 B1-30]" you get two chains.
    This is problematic because AF2 than folds this incorrectly as though D and B were also connected.
    To solve this chain IDs are changed using rechain.py.
    The script finds chainbreaks according to physical distance between CA atoms.
    """
    log.info("Running rechainRFdiffPDBs")
    rf_pdbs = glob.glob(os.path.join(cfg.rfdiff_out_path, "*.pdb"))
    for pdb in rf_pdbs:
        run_and_log(
            f'{cfg.pymol_python_path} {scripts_folder / "rechain.py"} "{pdb}" "{pdb}" --chain_break_cutoff_A {cfg.chain_break_cutoff_A}',
            cfg=cfg,
        )

        """
        When using potentials sometimes there were very disordered structures and a lot of chainbreaks
        This broke things down the line
        I am writing a fix for a symmetry case, do we need it for other things?
        If there are more chains then intended in symmetry, makes a folder for wrong oligomers and moves them there
        They're not processed further
        """
        if cfg.inference.symmetry != None:
            print(cfg.inference.symmetry)
            if "tetrahedral" in cfg.inference.symmetry:
                n_chains_intended = 4
            elif "icosahedral" in cfg.inference.symmetry:
                n_chains_intended = 20
            elif "octahedral" in cfg.inference.symmetry:
                n_chains_intended = 8
            else:
                n_chains_intended = int(cfg.inference.symmetry[1:])

            print(n_chains_intended)

            n_chains_actual = len(
                [i for i in PDBParser().get_structure("rechained", pdb).get_chains()]
            )

            if n_chains_actual != n_chains_intended:
                log.info(
                    f'intended N chains: {n_chains_intended}, actual N: {n_chains_actual}; moving structure to "disordered" folder'
                )
                disordered_dir = os.path.join(cfg.rfdiff_out_path, "disordered")
                print(disordered_dir)
                if not os.path.exists(disordered_dir):
                    os.makedirs(disordered_dir)
                shutil.move(pdb, disordered_dir)

    log.info("After rechaining")


def load_plugins(filter_configs: List[dict]):
    """
    Loads plugin scripts from a list of filter configs.
    Each config should include:
        - 'filter_name': a label
        - 'filter_script': path to .py file
        - 'delete_failed': whether to delete failed files
        - any kwargs specific to the plugin

    Returns:
        A list of (function, config) tuples, one per plugin.
    """
    plugins = []

    for config in filter_configs:
        script_path = (Path(__file__).parent / config["filter_script"]).resolve()
        log.debug(f"Resolved plugin path: {script_path}")
        if not script_path.exists():
            raise FileNotFoundError(f"Plugin script not found: {script_path}")

        # Dynamically load the module from the given path
        spec = importlib.util.spec_from_file_location(
            config["filter_name"], str(script_path)
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        # Make sure the plugin defines a filter_backbone function
        if not hasattr(module, "filter_backbone"):
            raise AttributeError(
                f"Plugin {config['filter_name']} does not define 'filter_backbone(pdb, kwargs)'"
            )

        # Store the function and the original config
        plugins.append((module.filter_backbone, config))
        log.info(f"Loaded plugin: {config['filter_name']} from {script_path}")

    return plugins


def filter_backbones_after_rfdiff(rfdiff_dir: Path, plugins):
    """
    Applies all filters defined in filter_configs to each PDB in the rfdiff output directory.

    Parameters:
        - rfdiff_dir: Path to the directory containing RFDiffusion output PDBs
        - filter_configs: List of plugin configs from the YAML config

    For each PDB:
        - If all plugins return True: keep the PDB
        - If any plugin returns False: either delete it or move it to 'failed_filters', based on config
    """
    rfdiff_dir = Path(rfdiff_dir)
    failed_dir = rfdiff_dir / "failed_filters"
    failed_dir.mkdir(exist_ok=True)  # Create a directory for failed files

    # plugins = load_plugins(filter_configs)  # Load the plugins

    pdbs = list(rfdiff_dir.glob("*.pdb"))  # Get all PDB files in the output dir
    if not pdbs:
        log.warning(f"No PDB files found in {rfdiff_dir}")
        return

    log.info(f"Applying filters to {len(pdbs)} PDB files...")

    for pdb_file in pdbs:
        try:
            passed_all = True
            trb_file = pdb_file.with_suffix(".trb")
            for filter_fn, config in plugins:
                # Only pass user-defined kwargs to plugin function
                kwargs = {
                    k: v
                    for k, v in config.items()
                    if k not in ["filter_name", "filter_script", "delete_failed"]
                }

                # Call the plugin’s filter_backbone function
                result = filter_fn(pdb_file, kwargs)

                if not result:
                    passed_all = False
                    if config.get("delete_failed", False):
                        pdb_file.unlink()  # delete file
                        if trb_file.exists():
                            trb_file.unlink()  # delete corresponding .trb file
                        log.info(
                            f"Deleted {pdb_file.name} and {trb_file.name} (failed '{config['filter_name']}')"
                        )
                    else:
                        dest = failed_dir / pdb_file.name
                        shutil.move(str(pdb_file), dest)  # move file
                        if trb_file.exists():
                            trb_file_dest = failed_dir / trb_file.name
                            shutil.move(
                                str(trb_file), trb_file_dest
                            )  # move corresponding .trb file
                        log.info(
                            f"Moved {pdb_file.name} and {trb_file.name} to {dest} (failed '{config['filter_name']}')"
                        )
                    break  # stop checking this PDB if it already failed one filter

            if passed_all:
                log.info(f"Passed: {pdb_file.name} — file kept in {rfdiff_dir}")

        except Exception as e:
            log.error(f"Error processing {pdb_file.name}: {e}")


def filter_backbone(backbone_file: Path, filter_configs: List[dict], plugins):
    """
    Applies all filters defined in filter_configs to each PDB in the rfdiff output directory.

    Parameters:
        - rfdiff_dir: Path to the directory containing RFDiffusion output PDBs
        - filter_configs: List of plugin configs from the YAML config

    For each PDB:
        - If all plugins return True: keep the PDB
        - If any plugin returns False: either delete it or move it to 'failed_filters', based on config
    """
    backbone_file = Path(backbone_file)
    rfdiff_dir = backbone_file.parent
    failed_dir = rfdiff_dir / "failed_filters"
    failed_dir.mkdir(exist_ok=True)  # Create a directory for failed files

    # plugins = load_plugins(filter_configs)  # Load the plugins

    log.info(f"Applying filters to {backbone_file}")

    try:
        passed_all = True
        trb_file = backbone_file.with_suffix(".trb")
        for filter_fn, config in plugins:
            # Only pass user-defined kwargs to plugin function
            kwargs = {
                k: v
                for k, v in config.items()
                if k not in ["filter_name", "filter_script", "delete_failed"]
            }

            # Call the plugin’s filter_backbone function
            result = filter_fn(backbone_file, kwargs)

            if not result:
                passed_all = False
                if config.get("delete_failed", False):
                    backbone_file.unlink()  # delete file
                    if trb_file.exists():
                        trb_file.unlink()  # delete corresponding .trb file
                    log.info(
                        f"Deleted {backbone_file.name} and {trb_file.name} (failed '{config['filter_name']}')"
                    )
                else:
                    dest = failed_dir / backbone_file.name
                    shutil.move(str(backbone_file), dest)  # move file
                    if trb_file.exists():
                        trb_file_dest = failed_dir / trb_file.name
                        shutil.move(
                            str(trb_file), trb_file_dest
                        )  # move corresponding .trb file
                    log.info(
                        f"Moved {backbone_file.name} and {trb_file.name} to {dest} (failed '{config['filter_name']}')"
                    )
                break  # stop checking this PDB if it already failed one filter

        if passed_all:
            log.info(f"Passed: {backbone_file.name} — file kept in {rfdiff_dir}")

    except Exception as e:
        log.error(f"Error processing {backbone_file.name}: {e}")


def parse_additional_args(cfg, group):
    dodatniArgumenti = ""
    for k, v in (cfg.get(group, {}) or {}).items():  # or to allow for empty groups
        dodatniArgumenti += f" {k} {v}"
    return dodatniArgumenti


error_messages = ["Testing if we can restart the prosculptApp(cfg)"]
import random


def throw(id, cycle=0, cfg="KEKEKE"):
    """
    For testing restartability. Call it at different points in the code to see if the restart works.
    @param id: int, ID of the crash you want to perform
    @param cycle: int, in which AF-MPNN cycle should it crash
    Message is derived from error_messages array.
    Alternatively, use probability-based crash test, which also re-runs the hard-coded command.

    To make it crash, add +throw=id +crash_at_cycle=cycle to cmd, before -cd and -cn
    """
    log.info(
        f"Wanna crash? {id} ?= {crash_at_error} in cycle {cycle} ?= {crash_at_cycle}"
    )
    # if id == crash_at_error and cycle == crash_at_cycle:
    if random.random() > 8 / 9:
        msg = error_messages[id]
        log.critical(f"Forced crash at {id}: {msg} in cycle {cycle}")
        log.info("Re-running the command: ")
        run_and_log(
            """python slurm_runner.py 1 MP_inpaintseq_CONT  output_dir="Examples/Examples_out/inpaintseq_MP_CONT" +throw=-16 +crash_at_cycle=0 -cd Examples -cn multipass_inpaintseq""",
            cfg=cfg,
        )
        raise Exception(f"Forced crash at {id}: {msg} in cycle {cycle}")


def save_checkpoint(folder, piece, value):
    with open(os.path.join(folder, f"checkpoint_{piece}.txt"), "w") as f:
        f.write(str(value))
        f.flush()
    log.info(f"Saving checkpoint {piece} = {value}")


def get_checkpoint(folder, piece, default=0):
    # check if file exists
    if os.path.exists(os.path.join(folder, f"checkpoint_{piece}.txt")):
        with open(os.path.join(folder, f"checkpoint_{piece}.txt"), "r") as f:
            value = int(f.read() or 0)
            log.info(f"Reading checkpoint {piece}: {value}")
            return value  # suppose only numbers (as strings) are there (user should not write to this file). In worst case, file is empty.
    else:
        log.info(f"Checkpoint {piece} doesn't exist. Returning default value {default}")
        return default


def plugin_filters(cfg):
    # Extract the filter config
    filter_config = cfg.get("rfdiff_backbone_filters", {})

    # Set the RFDiffusion output directory
    rfdiff_output_dir = cfg.rfdiff_out_dir  # Contains *.pdb files

    rfdiff_output_dir = Path(rfdiff_output_dir)
    plugins = load_plugins(filter_config)
    # Log the plugin configuration once
    log.info("Plugin configuration settings:")
    log.info(pprint.pformat(filter_config))
    # Run the plugin system
    filter_backbones_after_rfdiff(rfdiff_output_dir, plugins)

    # Check how many .pdb files remain
    remaining_pdbs = list(rfdiff_output_dir.glob("*.pdb"))
    if not remaining_pdbs:
        log.error(
            "All RFDiffusion backbones failed the plugin filters. No structures remain."
        )
        raise Exception(
            "All RFDiffusion backbones failed plugin filters — aborting pipeline."
        )


def do_cycling(cfg):
    """
    Runs a loop of (ProteinMPNN -> AF2 ->) `af2_mpnn_cycles` number of times.
    """
    log.info("Running do_cycling")
    start_cycle = get_checkpoint(cfg.output_dir, "cycle", 0)
    content_status = get_checkpoint(
        cfg.output_dir, "content_status", 0
    )  # 0 ... fresh run; 1 ... can clear and copy; 2 ... can copy; 3 ... corrupted (partial new files between MPNN and AF2); 4 ... AF2 is running, everything before AF2 should be skipped.
    for cycle in range(start_cycle, cfg.af2_mpnn_cycles):
        print("cycleeeeee", cycle)

        trb_paths = None
        input_mpnn = (
            cfg.rfdiff_out_dir
        )  # First cycle has this, other cycles overwrite this var

        if not cycle == 0:  # All other cycles get starting PDBs from AF2
            print(f"Cycle is not 0: {cycle}")
            cycle_directory = os.path.join(cfg.output_dir, "2_1_cycle_directory")

            if content_status == 1:
                if os.path.exists(cycle_directory):
                    shutil.rmtree(cycle_directory)  # We should not remove it on restart
                os.makedirs(cycle_directory, exist_ok=True)
                save_checkpoint(
                    cfg.output_dir, "content_status", 2
                )  ## AF folder is ok to move its files to cycle_folder and cycle_folder is empty (ready for new files).
                content_status = 2
            print(f"Nada: {cycle}")

            if content_status == 2:
                af2_model_subdicts = glob.glob(os.path.join(cfg.af2_out_dir, "*"))
                log.info(f"AF2 model subdicts: {af2_model_subdicts}")

                # Access all af2 models and put them in one intermediate directory to get the same ored as in the 1st cycle (all pdbs in one directory)
                for model_subdict in af2_model_subdicts:
                    log.info(f"Model subdict: {model_subdict}")
                    af2_pdbs = sorted(glob.glob(os.path.join(model_subdict, "T*.pdb")))
                    log.info(f"AF2 pdbs: {af2_pdbs}")

                    rf_model_num = prosculpt.get_token_value(
                        os.path.basename(model_subdict), "model_", "(\\d+)"
                    )
                    # Count number of files in cycle_directory which contain f"rf_{rf_model_num}"
                    number_of_copied_models = len(
                        [
                            f
                            for f in os.listdir(cycle_directory)
                            if f"rf_{rf_model_num}" in f
                        ]
                    )
                    log.info(
                        f"Number of already copied/present models: {number_of_copied_models}"
                    )
                    # last_iteration = max([int(prosculpt.get_token_value(f, "itr_", "(\\d+)")) for f in os.listdir(cycle_directory) if f"rf_{rf_model_num}" in f] or [0])

                    for i, af2_pdb in enumerate(
                        af2_pdbs, start=number_of_copied_models
                    ):
                        af_model_num = prosculpt.get_token_value(
                            os.path.basename(af2_pdb), "model_", "(\\d+)"
                        )
                        # if str(af_model_num) in args.af2_models:
                        # Rename pdbs to keep the model_num traceability with orginal rfdiff structure and enable filtering which models for next cycle

                        log.info(
                            f"i: {i}, af2_pdb: {af2_pdb}, af_model_num: {af_model_num}, rf_model_num: {rf_model_num}"
                        )
                        shutil.move(
                            af2_pdb,
                            os.path.join(
                                cycle_directory,
                                f"rf_{rf_model_num}__model_{af_model_num}__cycle_{cycle}__itr_{i}__.pdb",
                            ),
                        )
                        # rf_ --> rfdiffusion structure number (in rfdiff outou dir)
                        # model_ -> af2 model num, used for filtering which to cycle (preference for model 4)
                        # itr_ -> to differentiate models in later cycles (5 pdbs for model 4 from rf 0 for example)
                        # is it maybe possible to filter best ranked by af2 from the itr numbers?
                save_checkpoint(
                    cfg.output_dir, "content_status", 3
                )  ## ProteinMPNN is quick and can be rerun until AF2 starts running
                content_status = 3

            input_mpnn = cycle_directory

            if content_status != 4:
                shutil.rmtree(
                    cfg.mpnn_out_dir
                )  # Remove MPNN dir so you can create new sequences
                os.makedirs(cfg.mpnn_out_dir, exist_ok=True)  # Create it again
            trb_paths = os.path.join(cfg.rfdiff_out_dir, "*.trb")
            print("trb_path is: ", trb_paths)
        # endif

        mpnn_run_command_only_path = (cfg.mpnn_run_command).rsplit("/", 1)[0]

        if content_status != 4:
            # All cycles run the same commands
            run_and_log(
                f"{mpnn_run_command_only_path}/helper_scripts/parse_multiple_chains.py \
                --input_path={input_mpnn} \
                --output_path={cfg.path_for_parsed_chains}",
                cfg=cfg,
            )

            if cfg.chains_to_design:
                run_and_log(
                    f"{mpnn_run_command_only_path}/helper_scripts/assign_fixed_chains.py \
                        --input_path={cfg.path_for_parsed_chains} \
                        --output_path={cfg.path_for_assigned_chains} \
                        --chain_list='{cfg.chains_to_design}'",
                    cfg=cfg,
                )

            fixed_pos_path = prosculpt.process_pdb_files(
                input_mpnn, cfg.mpnn_out_dir, cfg, trb_paths, cycle=cycle
            )  # trb_paths is optional (default: None) and only used in non-first cycles
            # trb_paths is atm not used in process_pdb_files anyway -- a different approach is used (file.pdb -> withExtension .trb), which ensures the PDB and TRB files match.
            # but this function should still be called because it creates (probably) important .json files

            if cfg.inference.symmetry != None:
                run_and_log(
                    f"{mpnn_run_command_only_path}/helper_scripts/make_tied_positions_dict.py \
                    --input_path={cfg.path_for_parsed_chains} \
                    --output_path={cfg.path_for_tied_positions} \
                    --homooligomer 1",
                    cfg=cfg,
                )
                log.info(f"running symmetry")

            # _____________ RUN ProteinMPNN_____________
            # At first cycle, use num_seq_per_target from config. In subsequent cycles, set it to 1.
            proteinMPNN_cmd_str = f'{cfg.mpnn_run_command} \
                --jsonl_path {cfg.path_for_parsed_chains} \
                --fixed_positions_jsonl {cfg.path_for_fixed_positions} \
                {"--tied_positions_jsonl "+cfg.path_for_tied_positions if cfg.inference.symmetry!=None else ""} \
                --chain_id_jsonl {cfg.path_for_assigned_chains} \
                --out_folder {cfg.mpnn_out_dir} \
                --num_seq_per_target {cfg.num_seq_per_target_mpnn if cycle == 0 else 1} \
                --sampling_temp {cfg.sampling_temp} \
                --backbone_noise {cfg.backbone_noise} \
                --use_soluble_model  \
                --omit_AAs {cfg.omit_AAs} \
                {parse_additional_args(cfg, "pass_to_mpnn")} \
                --batch_size 1'

            run_and_log(proteinMPNN_cmd_str, cfg=cfg)

            log.info("Preparing to empty af2 directory.")

            # af2 directory must be empty
            shutil.rmtree(cfg.af2_out_dir)
            os.makedirs(cfg.af2_out_dir, exist_ok=True)

            # ________________ RUN AF2______________
            fasta_files = sorted(
                glob.glob(os.path.join(cfg.fasta_dir, "*.fa"))
            )  # glob is not sorted by default
            print(fasta_files)

            monomers_fasta_dir = os.path.join(cfg.fasta_dir, "monomers")
            if not os.path.exists(monomers_fasta_dir):
                os.makedirs(monomers_fasta_dir)

            # if symmetry - make fasta file with monomer sequence only
            if cfg.inference.symmetry != None or cfg.get(
                "model_monomer", False
            ):  # cfg symmetry is constant, not based on fasta_file
                for fasta_file in fasta_files:
                    sequences = []
                    for record in SeqIO.parse(fasta_file, "fasta"):
                        record.seq = record.seq[: record.seq.find("/")]
                        idd = record.id
                        record.id = "monomer_" + idd
                        descr = record.description
                        record.description = "monomer_" + descr
                        sequences.append(record)
                    SeqIO.write(
                        sequences,
                        os.path.join(
                            os.path.dirname(fasta_file),
                            "monomers",
                            f"{os.path.basename(fasta_file)[:-3]}_monomer.fa",
                        ),
                        "fasta",
                    )  # It must be written to the correct subfolder already to prevent duplication on job restart (CHECKPOINTING): _monomer_monomer.fa
                    print(
                        "File written to /monomers subfolder. " + fasta_file
                    )  # just for DEBUG
                    ## CHECKPOINTING: Right now, protMPNN is rerun after restart, which outputs different files into mpnn_out. The new files overwrite the old ones in /monomers, so it is always fresh.
                    ## Unles AF2 already started running -- then it is skipped.

        # content_status = 4; everything up to here can be skipped
        save_checkpoint(cfg.output_dir, "content_status", 4)
        content_status = 4

        fasta_files = glob.glob(os.path.join(cfg.fasta_dir, "*.fa"))
        fasta_files += glob.glob(os.path.join(cfg.fasta_dir, "monomers", "*.fa"))
        fasta_files = sorted(fasta_files)

        print(fasta_files)

        for (
            fasta_file
        ) in (
            fasta_files
        ):  # Number of fasta files corresponds to number of rfdiff models
            print("RUNNING FASTA ", fasta_file)
            if cycle == 0:
                rf_model_num = prosculpt.get_token_value(
                    os.path.basename(fasta_file), "_", "(\\d+)"
                )  # get 0 from _0.fa using reg exp
            else:
                rf_model_num = prosculpt.get_token_value(
                    os.path.basename(fasta_file), "rf_", "(\\d+)"
                )  # get 0 from rf_0__model_1__cycle_2__itr_0__.pdb
                af2_model_num = prosculpt.get_token_value(
                    os.path.basename(fasta_file), "model_", "(\\d+)"
                )  # get 1 from rf_0__model_1__cycle_2__itr_0__.pdb
            model_dir = os.path.join(
                cfg.af2_out_dir, f"model_{rf_model_num}"
            )  # create name for af2 directory name: model_0

            if "monomer" in os.path.basename(
                fasta_file
            ):  # if we are doing symmetry - make an extra directory for modeling monomers with af2
                model_dir = os.path.join(model_dir, "monomers")

            prosculpt.change_sequence_in_fasta(cfg.rfdiff_pdb, fasta_file)
            if not os.path.exists(model_dir):
                os.makedirs(model_dir)

            model_order = str(cfg.model_order).split(
                ","
            )  # If only one number is passed, Hydra converts it to int
            num_models = len(model_order)  # model order "1,2,3" -> num of models = 3

            if cfg.get("colabfold_preparation_script", False):
                preparation_command = f"source {cfg.colabfold_preparation_script} && "  # if the preparation script is specified, include it
            else:
                preparation_command = ""

            if cfg.get(
                "use_a3m", False
            ):  # If we're using a custom a3m, generate it for each sequence in the fasta
                if cfg.prediction_model == "Colabfold":
                    print("Generating custom msa files")
                    input_a3m_files = []
                    with open(fasta_file) as fasta_f:
                        a3m_filename = ""
                        for line in fasta_f:
                            if line[0] == ">":
                                a3m_filename = (
                                    "".join(
                                        [
                                            c if (c.isalnum() or c == ".") else "_"
                                            for c in line[1:-1]
                                        ]
                                    )
                                ) + ".a3m"
                            else:
                                mpnn_seq = line
                                trb_file = os.path.join(
                                    cfg.rfdiff_out_dir, "_" + str(rf_model_num) + ".trb"
                                )  #
                                print(a3m_filename)
                                custom_a3m_path = os.path.join(model_dir, a3m_filename)
                                prosculpt.make_alignment_file(
                                    cfg,
                                    trb_file,
                                    cfg.pdb_path,
                                    mpnn_seq,
                                    cfg.a3m_dir,
                                    custom_a3m_path,
                                )
                                input_a3m_files.append(custom_a3m_path)

                    for a3m_file in input_a3m_files:
                        if cycle == 0:  # have to run af2 differently in first cycle
                            run_and_log(
                                f'{preparation_command}{cfg.colabfold_run_command} \
                                            --model-type alphafold2_multimer_v3 \
                                            --msa-mode mmseqs2_uniref_env \
                                            {a3m_file} {model_dir} \
                                            --model-order {cfg.model_order} \
                                            {parse_additional_args(cfg, "pass_to_af")} \
                                            --num-models {num_models}',
                                cfg=cfg,
                            )  # changed from single_sequence
                        else:
                            for model_number in model_order:
                                if (
                                    af2_model_num == model_number
                                ):  # From the af2 model 4 you want only model 4 not also 2 and for 2 only 2 not 4 (--model_order "2,4")
                                    num_models = 1
                                    run_and_log(
                                        f'{preparation_command}{cfg.colabfold_run_command} \
                                            --model-type alphafold2_multimer_v3 \
                                            --msa-mode mmseqs2_uniref_env \
                                            {a3m_file} {model_dir} \
                                            --model-order {model_number} \
                                            {parse_additional_args(cfg, "pass_to_af")} \
                                            --num-models {num_models}',
                                        cfg=cfg,
                                    )  # Changed from single_sequence
                elif cfg.prediction_model == "Boltz2":  # If using Boltz
                    print("Generating custom msa files for Boltz")

                    yaml_dir = os.path.join(model_dir, "yaml_inputs")
                    os.makedirs(yaml_dir, exist_ok=True)
                    alignment_inputs_dir = os.path.join(model_dir, "alignment_inputs")
                    os.makedirs(alignment_inputs_dir, exist_ok=True)
                    input_yaml_files = []
                    with open(fasta_file) as fasta_f:
                        a3m_filename = ""
                        for line in fasta_f:
                            if line[0] == ">":
                                sequence_id = "".join(
                                    [
                                        c if (c.isalnum() or c == ".") else "_"
                                        for c in line[1:-1]
                                    ]
                                )
                            else:
                                mpnn_seq = line
                                trb_file = os.path.join(
                                    cfg.rfdiff_out_dir, "_" + str(rf_model_num) + ".trb"
                                )  #
                                print(a3m_filename)

                                prosculpt.make_alignment_file_boltz(
                                    sequence_id,
                                    mpnn_seq,
                                    cfg.a3m_dir,
                                    alignment_inputs_dir,
                                )
                                custom_yaml_path = prosculpt.make_boltz_input_yaml(
                                    cfg,
                                    sequence_id,
                                    mpnn_seq,
                                    yaml_dir,
                                    alignment_inputs_dir,
                                )
                                input_yaml_files.append(custom_yaml_path)

                    run_and_log(
                        f"boltz predict {yaml_dir}/ --out_dir {model_dir} --output_format pdb {cfg.get('boltz_options', '')}",
                        cfg=cfg,
                    )
                else:
                    log.error(f"Unsupported prediction model: {cfg.prediction_model}")
            else:  # this is the normal mode of operations. Single sequence.
                if cfg.prediction_model == "Colabfold":
                    if cycle == 0:  # have to run af2 differently in first cycle
                        run_and_log(
                            f'{preparation_command}{cfg.colabfold_run_command} \
                                        --model-type alphafold2_multimer_v3 \
                                        --msa-mode single_sequence \
                                        {fasta_file} {model_dir} \
                                        --model-order {cfg.model_order} \
                                        {parse_additional_args(cfg, "pass_to_af")} \
                                        --num-models {num_models}',
                            cfg=cfg,
                        )
                    else:
                        for model_number in model_order:
                            if (
                                af2_model_num == model_number
                            ):  # From the af2 model 4 you want only model 4 not also 2 and for 2 only 2 not 4 (--model_order "2,4")
                                num_models = 1
                                run_and_log(
                                    f'{preparation_command}{cfg.colabfold_run_command} \
                                        --model-type alphafold2_multimer_v3 \
                                        --msa-mode single_sequence \
                                        {fasta_file} {model_dir} \
                                        --model-order {model_number} \
                                        {parse_additional_args(cfg, "pass_to_af")} \
                                        --num-models {num_models}',
                                    cfg=cfg,
                                )
                elif cfg.prediction_model == "Boltz2":  # If using Boltz
                    yaml_dir = os.path.join(model_dir, "yaml_inputs")
                    os.makedirs(yaml_dir, exist_ok=True)

                    with open(fasta_file) as fasta_f:
                        for line in fasta_f:
                            if line[0] == ">":
                                sequence_id = "".join(
                                    [
                                        c if (c.isalnum() or c == ".") else "_"
                                        for c in line[1:-1]
                                    ]
                                )
                            else:
                                mpnn_seq = line

                                custom_yaml_path = prosculpt.make_boltz_input_yaml(
                                    cfg,
                                    sequence_id,
                                    mpnn_seq,
                                    yaml_dir,
                                    None,
                                )
                                # input_yaml_files.append(custom_yaml_path)

                    run_and_log(
                        f"boltz predict {yaml_dir}/ --out_dir {model_dir} --output_format pdb {cfg.get('boltz_options', '')}",
                        cfg=cfg,
                    )
                else:
                    log.error(f"Unsupported prediction model: {cfg.prediction_model}")
        save_checkpoint(
            cfg.output_dir, "content_status", 1
        )  ## Content is ok to be copied to cycle_dir
        save_checkpoint(
            cfg.output_dir, "cycle", cycle + 1
        )  # Restart should start with next cycle.
        content_status = 1

        # msa single sequence makes sense for designed proteins


def final_operations(cfg):
    log.info("Final operations")
    json_directories = glob.glob(os.path.join(cfg.af2_out_dir, "*"))

    print(
        "do we already have csv files? If we need to re-run stats, we have to delete them: ",
        os.listdir(cfg.output_dir),
    )
    if "output.csv" in os.listdir(cfg.output_dir):
        os.remove(cfg.output_dir + "/output.csv")
    if "rosetta_scores.csv" in os.listdir(cfg.output_dir):
        os.remove(cfg.output_dir + "/rosetta_scores.csv")
    if "final_output.csv" in os.listdir(cfg.output_dir):
        os.remove(cfg.output_dir + "/final_output.csv")

    for model_i in json_directories:  # for model_i in [model_0, model_1, model_2 ,...]

        trb_num = prosculpt.get_token_value(
            os.path.basename(model_i), "model_", "(\\d+)"
        )  # get 0 from model_0 using reg exp
        if cfg.prediction_model == "Colabfold":
            if "pdb_path" in cfg:
                prosculpt.rename_pdb_create_csv_colabfold(
                    cfg,
                    cfg.output_dir,
                    cfg.rfdiff_out_dir,
                    trb_num,
                    model_i,
                    cfg.pdb_path,
                    cfg.inference.symmetry,
                    model_monomer=cfg.get("model_monomer", False),
                )
            else:
                prosculpt.rename_pdb_create_csv_colabfold(
                    cfg,
                    cfg.output_dir,
                    cfg.rfdiff_out_dir,
                    trb_num,
                    model_i,
                    control_structure_path=None,
                    symmetry=cfg.inference.symmetry,
                    model_monomer=cfg.get("model_monomer", False),
                )
        elif cfg.prediction_model == "Boltz2":  # If using Boltz
            if "pdb_path" in cfg:
                prosculpt.rename_pdb_create_csv_boltz(
                    cfg,
                    cfg.output_dir,
                    cfg.rfdiff_out_dir,
                    trb_num,
                    model_i,
                    cfg.pdb_path,
                    cfg.inference.symmetry,
                    model_monomer=cfg.get("model_monomer", False),
                )
            else:
                prosculpt.rename_pdb_create_csv_boltz(
                    cfg,
                    cfg.output_dir,
                    cfg.rfdiff_out_dir,
                    trb_num,
                    model_i,
                    control_structure_path=None,
                    symmetry=cfg.inference.symmetry,
                    model_monomer=cfg.get("model_monomer", False),
                )
        else:
            log.error(f"Unsupported prediction model: {cfg.prediction_model}")
    csv_path = os.path.join(
        cfg.output_dir, "output.csv"
    )  # constructed path 'output.csv defined in rename_pdb_create_csv function
    run_and_log(
        f'{cfg.prosculpt_python_path} {scripts_folder / "scoring_script.py"} {csv_path}',
        cfg=cfg,
    )

    rosetta_scores_path = os.path.join(
        cfg.output_dir, "rosetta_scores.csv"
    )  #'rosetta_scores.csv defined in scoring_rg_... script
    prosculpt.merge_csv(
        cfg.output_dir, csv_path, rosetta_scores_path
    )  # , cfg.get("output_best", True), cfg.get("rmsd_threshold",3),cfg.get("plddt_threshold",90) used for filtering. not needed now.

    os.remove(csv_path)
    os.remove(rosetta_scores_path)


def pass_config_to_rfdiff(cfg):
    """
    This function creates a new .yaml file with "inference", "potentials" etc. read from run.yaml. Also add defaults: - base \n - _self_
    Then, when calling RfDiff, add -cd [folder_with_this_created_config] -cn [name_of_created_config_file]
    """

    # Keep in mind that if passing through cmd, these parameters need to have ++ in front of them.
    # TODO: should everything go through the command line? Makes it a bit cleaner.
    def keep_only_groups(pair):
        groupsToPassToRfDiff = (
            cfg.pass_to_rfdiff
        )  # groups that should be passed to RfDiff
        key, value = pair
        if key not in groupsToPassToRfDiff:
            return False
        else:
            return True

    new_cfg = dict(filter(keep_only_groups, cfg.items()))
    new_cfg["defaults"] = [
        "base",
        "_self_",
    ]  # TODO: is the base always the right thing to include? For symm it should be symmetry.
    if "pdb_path" in cfg:
        with open_dict(new_cfg["inference"]):
            new_cfg["inference"]["input_pdb"] = cfg.pdb_path

    log.info(f"Saving this config for RfDiff: {OmegaConf.to_yaml(new_cfg)}")
    # save to file
    with open(os.path.join(cfg.output_dir, "prosculpt2rfdiff.yaml"), "w") as f:
        f.write(OmegaConf.to_yaml(new_cfg))

    pass


import time

TIMEMEASURES = {}
TIMECALC = time.time()
PREVIOUS_MESSAGE = "Before app start"


def dtimelog(message, final=False):
    global TIMECALC
    global PREVIOUS_MESSAGE
    dt = round(time.time() - TIMECALC, 1)
    log.info(f"* * * {PREVIOUS_MESSAGE} lasted {dt} s. Running {message} * * *")
    TIMEMEASURES[PREVIOUS_MESSAGE] = dt
    PREVIOUS_MESSAGE = message
    TIMECALC = time.time()

    if final:
        # Print all measures one per line
        for k, v in TIMEMEASURES.items():
            log.info(f"{k} lasted {v} s.")


crash_at_error = 0  # Added through command line. Where should the test crash?
crash_at_cycle = 0  # Added through command line. In which cycle should the test crash?


@hydra.main(version_base=None, config_path="config", config_name="run")
def prosculptApp(cfg: DictConfig) -> None:

    log.info("Hydrić")
    log.info("The following configuration was passed: \n" + OmegaConf.to_yaml(cfg))

    # log.debug(f"Now in Hydra, cwd = {os.getcwd()}")

    dtimelog("general_config_prep")
    general_config_prep(cfg)
    config = HydraConfig.get()
    config_name = config.job.config_name
    log.info(f"config.runtime.config_sources: {config.runtime.config_sources}")
    config_path = [
        path["path"]
        for path in config.runtime.config_sources
        if path["schema"] == "file"
    ][
        -1
    ]  # If run without config (through cmd), [1] is out of range. Assume passed config file is always the last one.
    shutil.copy(
        os.path.join(config_path, config_name + ".yaml"),
        os.path.join(cfg.output_dir, f"input.yaml"),
    )

    global crash_at_error
    global crash_at_cycle
    crash_at_error = cfg.get("throw", -1)
    crash_at_cycle = cfg.get("crash_at_cycle", 0)

    if cfg.get("only_run_analysis", False):  # only_run_analysis
        print("***Skip everything, go to final operations***")
        dtimelog("only final_operations")
        final_operations(cfg)

    else:
        if not cfg.get("skipRfDiff", False):
            dtimelog("pass_config_to_rfdiff")
            pass_config_to_rfdiff(cfg)
            dtimelog("run_rfdiff")
            run_rfdiff(cfg)

            dtimelog("plugin_filters")
            plugin_filters(cfg)

            dtimelog("rechain_rfdiff_pdbs")
            rechain_rfdiff_pdbs(cfg)
        else:
            log.info("*** Skipping RfDiff ***")
            dtimelog("skipping RfDiff")
            # Copy input PDB to RfDiff_output_dir and rename it to follow the token scheme
            shutil.copy(cfg.pdb_path, os.path.join(cfg.rfdiff_out_dir, "_0.pdb"))

        dtimelog("do_cycling")
        do_cycling(cfg)
        dtimelog("final_operations")
        final_operations(cfg)
        dtimelog("Finished", True)


if __name__ == "__main__":
    print("File: ", __file__)
    a = pathlib.Path(__file__).resolve().parent
    log.info(f"Before running MYapp, cwd = {os.getcwd()}")
    prosculptApp()
