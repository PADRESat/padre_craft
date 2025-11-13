"""Provides functions to upload data to the time series database for display"""

from astropy.timeseries import TimeSeries
from padre_meddea.housekeeping.calibration import calibrate_hk_ts
from swxsoc.util.util import record_timeseries

import padre_craft.util.util as util


def record_housekeeping(hk_ts: TimeSeries, data_type):
    """Send the housekeeping time series to AWS."""
    if data_type == "meddea":
        hk_ts = util.convert_meddea_colnames(hk_ts)
        cal_hk_ts = calibrate_hk_ts(hk_ts)
        col_to_removes = ["CCSDS1", "CCSDS3", "checksum"]
        for this_col in col_to_removes:
            cal_hk_ts.remove_column(this_col)
    record_timeseries(hk_ts, data_type, "craft")
