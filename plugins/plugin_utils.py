import logging

logger = logging.getLogger("prosculpt.plugins")

_pyrosetta_initialized = False  # Global flag to avoid double init


def ensure_pyrosetta_initialized():
    """
    Safely initialize PyRosetta once globally.
    If already initialized, does nothing.
    If PyRosetta isn't installed, logs a warning.
    """
    global _pyrosetta_initialized

    if _pyrosetta_initialized:
        logger.info("PyRosetta already initialized.")
        return

    try:
        logger.info("Attempting to import pyrosetta...")
        from pyrosetta import init
        init(extra_options="-mute all")
        _pyrosetta_initialized = True
        logger.info("PyRosetta initialized.")
    except ImportError:
        logger.warning("PyRosetta not installed — skipping initialization.")
