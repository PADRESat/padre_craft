"""Provides functions to upload data to the time series database for display"""

from astropy.timeseries import TimeSeries
from swxsoc.util.util import record_timeseries

import padre_craft.util.util as util


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
    record_timeseries(my_ts, data_type, "craft")


def record_orbit(padre_orbit_ts: TimeSeries):
    """Send the orbit time series to AWS."""
    record_timeseries(padre_orbit_ts, "orbit", "craft")
