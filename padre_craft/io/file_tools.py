"""
Provides generic file readers.
"""

from pathlib import Path

import astropy.io.fits as fits
from astropy.io import ascii
from astropy.time import Time
from astropy.timeseries import TimeSeries

from padre_craft import LAUNCH_DATE
from padre_craft.util.util import filename_to_datatype

__all__ = ["read_file", "read_raw_file", "read_fits"]


def read_file(filename: Path):
    """
    Read a file.

    Parameters
    ----------
    filename: Path
        A file to read.

    Returns
    -------
    data

    Examples
    --------
    """
    # TODO: the following should use parse_science_filename
    this_path = Path(filename)
    match this_path.suffix.lower():
        case ".csv":  # raw binary file
            result = read_raw_file(this_path)
        case ".fits":  # level 0 or above
            result = read_fits(this_path)
        case _:
            raise ValueError(f"File extension {this_path.suffix} not recognized.")
    return result


def read_raw_file(file_path: Path) -> TimeSeries:
    """
    Read a raw (csv) data file and return a timeseries.
    Note that columns that cannot are not recognized as ints or floats are removed.

    Parameters
    ---------
    filename : Path
        A file to read

    Returns
    -------
    ts : TimeSeries
        The timeseries data
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)

    # Read in the CSV
    data_table = ascii.read(file_path, format="csv")

    # Make sure the time column exists
    time_column_name = "timestamp_ms"
    if time_column_name not in data_table.colnames:
        raise ValueError("No time column found, timestamp_ms")

    # Convert Unix milliseconds to Astropy Time
    time = Time(data_table[time_column_name] / 1000.0, format="unix")
    time.format = "isot"

    # Create the TimeSeries
    ts = TimeSeries(time=time, data=data_table)

    # Validate Times in the Data
    if any(ts.time < LAUNCH_DATE) > 0:
        raise ValueError("Found time before launch.")

    # Filter out non-numeric columns except time
    bad_col_names = []
    for this_col in ts.itercols():
        if not isinstance(this_col, Time):
            if this_col.dtype.kind not in ["i", "f"]:
                bad_col_names.append(this_col.name)
    if len(bad_col_names) > 0:
        ts.remove_columns(bad_col_names)

    # Record Metadata about the File
    ts.meta.update({"filename": file_path.name})
    ts.meta.update({"data_type": filename_to_datatype(file_path.name)})
    ts.sort()
    return ts


def read_fits(filename: Path):
    """
    Read a fits file of any level and return the appropriate data objects.
    """
    hdul = fits.open(filename)
    # header = hdul[0].header.copy()

    return hdul
