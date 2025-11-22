from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time
from astropy.timeseries import TimeSeries
from swxsoc.util import create_science_filename, parse_science_filename

from padre_craft import log

__all__ = [
    "filename_to_datatype",
    "create_craft_filename",
    "parse_science_filename",
    "convert_meddea_colnames",
]

TOKEN_TO_DATATYPE = {
    "CUBEADCS": "adcs",
    "EPS": "housekeeping",
    "GNSS": "gnss",
    "MEDDEA": "meddea",
    "SHIP": "sharp",
    "BP": "battery",
    "OBC_0": "obc",
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


def convert_meddea_colnames(ts: TimeSeries):
    """Given the column names from OBC MeDDEA standard to padre_meddea standard."""
    # translation betweem OBC MeDDEA housekeeping names to padre_meddea HK names

    OBC_TO_MEDDEA = {
        "FPTemp": "fp_temp",
        "DIBTemp": "dib_temp",
        "HVTemp": "hvps_temp",
        "HVVolts": "hvps_vsense",
        "HVCurrent": "hvps_csense",
        "Amps_1V5": "csense_15v",
        "Amps_3V3_D": "csense_33vd",
        "Amps_3V3_A": "csense_33va",
        "phRate": "hit_rate",
        "goodCmdCount": "good_cmd_cnt",
        "errorCount": "error_cnt",
        "heaterPWM": "heater_pwm_duty_cycle",
        "decimationRate": "decimation_rate",
        "sysError": "error_summary",
    }

    for obc_col, meddea_col in OBC_TO_MEDDEA.items():
        if obc_col in ts.colnames:
            ts.rename_column(obc_col, meddea_col)
    return ts


def remove_bad_data(ts: TimeSeries) -> TimeSeries:
    """Given a padre craft timeseries, remove all rows with bad data. Checks the following conditions
    * all data rows sum to zero.
    * remove rows with times that are in the future
    Returns
    -------
    filtered_ts
    """
    tbl = Table(ts)
    tbl.remove_column("time")
    if "timestamp_ms" in tbl.colnames:
        tbl.remove_column("timestamp_ms")
    row_sum = np.sum([tbl[col] for col in tbl.colnames], axis=0)
    good_times = ts.time <= Time.now()
    filtered_index = (row_sum != 0) * good_times
    ts_filtered = ts[filtered_index]
    if len(ts_filtered) != len(ts):
        log.warning(
            f"Found {len(ts) - len(ts_filtered)} rows of bad data and removed them."
        )
    return ts_filtered


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
