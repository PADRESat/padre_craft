"""Provides functions to upload data to the time series database for display"""

from astropy.timeseries import TimeSeries
from padre_meddea.housekeeping.calibration import calibrate_hk_ts
from swxsoc.util.util import record_timeseries

import padre_craft.util.util as util
from padre_craft.dirlist.dirlist import DirList


def record_housekeeping(hk_ts: TimeSeries, data_type: str) -> None:
    """
    Record housekeeping time series data to AWS TimeStream database.
    This function processes and stores housekeeping telemetry data, with special
    handling for MEDDEA data.

    Parameters
    ----------
    hk_ts : TimeSeries
        The housekeeping time series data to be recorded. The original time series
        is not modified; a copy is created for processing.
    data_type : str
        The type of data being recorded (e.g., "meddea"). This determines the
        processing pipeline applied to the data.
    """
    # Create a copy to avoid modifying the original
    my_ts = hk_ts.copy()

    # Special handling for MEDDEA data
    if data_type == "meddea":
        # Convert MEDDEA column names from OBC input to PADRE MEDDEA standard
        my_ts = util.convert_meddea_colnames(my_ts)
        # Remove unwanted columns
        col_to_removes = ["CCSDS1", "CCSDS3", "checksum", "timestamp_ms"]
        for this_col in col_to_removes:
            if this_col in hk_ts.colnames:
                my_ts.remove_column(this_col)

        # Fix Bad Data in the Time Series
        my_ts = util.remove_bad_data(my_ts)

        # Apply calibration from padre_meddea to housekeeping data
        my_ts = calibrate_hk_ts(my_ts)

    # Record the Time Series in the TimeStream Database
    record_timeseries(ts=my_ts, ts_name=data_type, instrument_name="craft")


def record_orbit(padre_orbit_ts: TimeSeries) -> None:
    """Send the orbit time series to AWS."""
    record_timeseries(padre_orbit_ts, "orbit", "craft")


def record_dirlist(this_dirlist: DirList) -> None:
    """
    Record directory listing summary data (file sizes and counts) to AWS.

    This function converts a DirList into summary time series for file sizes
    and file counts, then uploads those summaries to the time series database.

    Parameters
    ----------
    this_dirlist : DirList
        Directory listing to be summarized into file size and file count
        time series for upload.
    """
    # record the file sizes
    summary_ts = this_dirlist.to_summary_ts(type="size")
    record_timeseries(
        ts=summary_ts, ts_name="dirlist_file_size", instrument_name="craft"
    )
    # record the file counts
    summary_ts = this_dirlist.to_summary_ts(type="count")
    record_timeseries(
        ts=summary_ts, ts_name="dirlist_file_count", instrument_name="craft"
    )
