"""Provides functions to upload data to the time series database for display"""

from astropy.timeseries import TimeSeries
from padre_meddea.housekeeping.calibration import calibrate_hk_ts
from swxsoc.util.util import record_timeseries

import padre_craft.util.util as util


def record_housekeeping(hk_ts: TimeSeries, data_type):
    """Send the housekeeping time series to AWS."""
    my_ts = hk_ts.copy()
    if data_type == "meddea":
        my_ts = util.convert_meddea_colnames(my_ts)
        my_ts = calibrate_hk_ts(my_ts)
        col_to_removes = ["CCSDS1", "CCSDS3", "checksum"]
        for this_col in col_to_removes:
            if this_col in hk_ts.colnames:
                my_ts.remove_column(this_col)
    record_timeseries(my_ts, data_type, "craft")
