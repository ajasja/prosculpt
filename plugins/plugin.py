import importlib.util    
import logging            
import shutil              
from pathlib import Path   
from typing import List
import pprint


# Sets up global logging format and level.
logger = logging.getLogger(__name__)



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
        script_path = Path(config["filter_script"])
        if not script_path.exists():
            raise FileNotFoundError(f"Plugin script not found: {script_path}")

        # Dynamically load the module from the given path
        spec = importlib.util.spec_from_file_location(config["filter_name"], str(script_path))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        # Make sure the plugin defines a filter_backbone function
        if not hasattr(module, "filter_backbone"):
            raise AttributeError(f"Plugin {config['filter_name']} does not define 'filter_backbone(pdb, kwargs)'")

        # Store the function and the original config
        plugins.append((module.filter_backbone, config))
        logger.info(f"Loaded plugin: {config['filter_name']} from {script_path}")

    return plugins





def filter_backbones_after_rfdiff(rfdiff_dir: Path, filter_configs: List[dict]):
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

    plugins = load_plugins(filter_configs)  # Load the plugins

    # Log the plugin configuration once
    logger.info("Plugin configuration settings:")
    logger.info(pprint.pformat(filter_configs))

    pdbs = list(rfdiff_dir.glob("*.pdb"))  # Get all PDB files in the output dir
    if not pdbs:
        logger.warning(f"No PDB files found in {rfdiff_dir}")
        return

    logger.info(f"Applying filters to {len(pdbs)} PDB files...")

    for pdb_file in pdbs:
        try:
            passed_all = True

            for filter_fn, config in plugins:
                # Only pass user-defined kwargs to plugin function
                kwargs = {
                    k: v for k, v in config.items()
                    if k not in ["filter_name", "filter_script", "delete_failed"]
                }

                # Call the plugin’s filter_backbone function
                result = filter_fn(pdb_file, kwargs)

                if not result:
                    passed_all = False
                    if config.get("delete_failed", False):
                        pdb_file.unlink()  # delete file
                        logger.info(f"Deleted {pdb_file.name} (failed '{config['filter_name']}')")
                    else:
                        dest = failed_dir / pdb_file.name
                        shutil.move(str(pdb_file), dest)  # move file
                        logger.info(f"Moved {pdb_file.name} to {dest} (failed '{config['filter_name']}')")
                    break  # stop checking this PDB if it already failed one filter

            if passed_all:
                logger.info(f"Passed: {pdb_file.name} — file kept in {rfdiff_dir}")

        except Exception as e:
            logger.error(f"Error processing {pdb_file.name}: {e}")

