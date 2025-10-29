from pathlib import Path

from astropy.time import Time
from swxsoc.util import (
    create_science_filename,
    parse_science_filename,
)

from padre_craft import log

__all__ = ["filename_to_datatype", "create_craft_filename", "parse_science_filename"]

TOKEN_TO_DATATYPE = {
    "CUBEADCS": "adcs",
    "EPS": "housekeeping",
    "GNSS": "gnss",
    "MEDDEA": "meddea",
    "SHIP": "sharp",
    "BP": "battery",
}


def filename_to_datatype(filename: Path) -> str:
    """Given a filename path, return the data type"""
    if not isinstance(filename, Path):
        filename = Path(filename)
    token = filename.name.split("get_")[1].split("_Data")[0]
    for this_str, datatype in TOKEN_TO_DATATYPE.items():
        if this_str in token:
            return datatype
    log.warning(f"Could not determine data type for file {filename.name}")
    return token


def create_craft_filename(
    time: Time,
    level: str,
    descriptor: str,
    version: str,
    test: bool = False,
    overwrite: bool = False,
) -> str:
    """
    Generate the MEDDEA filename based on the provided parameters.

    Parameters
    ----------
    time : Time
        The time associated with the data.
    level : str
        The data level (e.g., "L1", "L2").
    descriptor : str
        The data descriptor (e.g., "SCI", "CAL").
    test : str
        The test identifier (e.g., "TEST1", "TEST2").
    overwrite : bool
        Whether to overwrite existing files.

    Returns
    -------
    str
        The generated MEDDEA filename.
    """
    # Filename Version X.Y.Z comes from two parts:
    #   1. Files Version Base: X.Y comes from the Software Version -> Data Version Mapping
    #   2. File Version Incrementor: Z starts at 0 and iterates for each new version based on what already exists in the filesystem.
    # version_base = "1.0"
    # version_increment = 0
    # version_str = f"{version_base}.{version_increment}"
    version_str = version

    # The Base Filename is used for searching to see if we need to increase our version increment.
    base_filename = create_science_filename(
        instrument="meddea",
        time=time,
        level=level,
        descriptor=descriptor,
        test=test,
        version=version_str,
    )
    base_filename = base_filename.replace("meddea", "craft")
    return base_filename
