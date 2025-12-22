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
    """
    Convert a filename to its corresponding data type descriptor.

    This function extracts the data type from a filename by parsing tokens between
    'get_' and '_Data' in the filename, then matches it against known data type
    mappings in TOKEN_TO_DATATYPE.

    Parameters
    ----------
    filename : Path or str
        The path to the file whose data type needs to be determined. If a string
        is provided, it will be converted to a Path object.

    Returns
    -------
    str
        The data type corresponding to the filename if found in TOKEN_TO_DATATYPE,
        otherwise returns the parsed token with a warning logged.

    Examples
    --------
    >>> filename_to_datatype(Path("padre_get_MEDDEA_HOUSE_KEEPING_Data_1762493454480_1762611866270.csv"))
    'meddea'

    Notes
    -----
    The function expects filenames to follow the pattern: *get_<descriptor>_Data*
    where <descriptor> contains the data type identifier that can be mapped to
    a known data type via TOKEN_TO_DATATYPE dictionary.

    Warnings
    --------
    If no matching data type is found in TOKEN_TO_DATATYPE, a warning is logged
    and the raw parsed token is returned.
    """
    if not isinstance(filename, Path):
        filename = Path(filename)
    # Parse out the "Descriptor" from the filename
    token = filename.name.split("get_")[1].split("_Data")[0]
    # Search for known "Descriptors" in the parsed token
    for this_str, datatype in TOKEN_TO_DATATYPE.items():
        if this_str in token:
            return datatype
    log.warning(f"Could not determine data type for file {filename.name}")
    return token


def convert_meddea_colnames(ts: TimeSeries) -> TimeSeries:
    """
    Convert MeDDEA column names from OBC standard to padre_meddea standard.
    This function renames columns in a TimeSeries object from the OBC (Onboard Computer)
    MeDDEA housekeeping naming convention to the PADRE MeDDEA housekeeping naming convention.
    Only columns that exist in the input TimeSeries will be renamed.

    Parameters
    ----------
    ts : TimeSeries
        A TimeSeries object containing MeDDEA data with OBC standard column names.
    Returns
    -------
    TimeSeries
        The same TimeSeries object with renamed columns following the padre_meddea standard.
        Columns not listed in the mapping dictionary remain unchanged.
    """
    # translation between OBC MeDDEA housekeeping names to padre_meddea HK names

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
    """
    Remove bad data from a PADRE craft TimeSeries by identifying and setting invalid rows to NaN.
    This function identifies rows with invalid data based on two conditions:
    1. All data columns (excluding time) sum to zero
    2. Timestamps are in the future (beyond current time)
    Rows meeting either condition are flagged as bad data, logged, and all their non-time
    column values are set to NaN. The time column is preserved.

    Parameters
    ----------
    ts : TimeSeries
        Input PADRE craft timeseries containing time and data columns to be filtered.

    Returns
    -------
    ts : TimeSeries
        The input timeseries with bad data rows set to NaN.
    """
    # Convert Astropy TimeSeries to Astropy Table for easier manipulation
    tbl = Table(ts)
    # Remove Time column
    tbl.remove_column("time")

    # Create a boolean filter for rows with all zero data or future timestamps
    row_sum = np.sum([tbl[col] for col in tbl.colnames], axis=0)
    good_times = ts.time <= Time.now()
    filtered_index = (row_sum != 0) * good_times

    # For values within the Filter, set all values to NaN
    # Log the bad data rows before setting to NaN
    bad_data_mask = ~filtered_index
    if np.any(bad_data_mask):
        bad_indices = np.where(bad_data_mask)[0]
        log.warning(f"Found {len(bad_indices)} rows of bad data. Setting to NaN.")

        # Set bad data rows to NaN
        for col in ts.colnames:
            if col != "time":  # Don't modify the time column
                # Convert TimeSeries Columns to float data types if necessary
                ts[col] = ts[col].astype(float)
                # Set Bad Data to NaN
                ts[col][bad_data_mask] = np.nan

    return ts


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
