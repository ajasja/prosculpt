from pathlib import Path
import logging

from pyrosetta import pose_from_file
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

logger = logging.getLogger(__name__)

def get_ss_xml():
    return """
<ROSETTASCRIPTS>

    <FILTERS>
        <SecondaryStructureCount
            name="ss_count"
            filter_helix_sheet="1"
            num_helix_sheet="3"
            return_total="1"
        />
    </FILTERS>

    <PROTOCOLS>
        <Add filter_name="ss_count" />
    </PROTOCOLS>

</ROSETTASCRIPTS>
"""

# Plugin entry point
def filter_backbone(pdb: str | Path, kwargs: dict):
    """
    Passes if the structure contains more than 2 DSSP-defined
    secondary structure elements (helices + beta strands)
    """

    logger.info(f"Running secondary structure element filter on: {pdb}")

    pose = pose_from_file(str(pdb))

    xml_obj = XmlObjects.create_from_string(get_ss_xml())
    ss_filter = xml_obj.get_filter("ss_count")

    passed = ss_filter.apply(pose)

    # Get the actual SS element count
    try:
        ss_count = ss_filter.report_sm(pose)
    except Exception:
        ss_count = "N/A"

    logger.info(
        f"Secondary structure elements: {ss_count} | Filter passed: {passed}"
    )

    return passed