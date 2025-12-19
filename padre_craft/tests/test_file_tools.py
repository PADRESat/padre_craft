import numpy as np
import pytest
from astropy.time import Time
from astropy.timeseries import TimeSeries

import padre_craft.io.file_tools as file_tools
from padre_craft import _test_files_directory

test_file_paths = _test_files_directory.glob("padre_*.csv")


@pytest.mark.parametrize("this_path", list(test_file_paths))
def test_file_read(this_path):
    """Test that all test files can be read"""
    ts = file_tools.read_file(this_path)
    assert isinstance(ts, TimeSeries)
    # check that there are no unexpected column data types
    for this_col in ts.itercols():
        if isinstance(this_col, Time):
            continue
        else:
            assert this_col.dtype.kind in ["i", "f"]
    if "MEDDEA" in this_path.name:
        assert len(ts) == 13
        assert len(np.unique(ts.time)) == 13
    else:
        assert len(ts) == 10
        assert len(np.unique(ts.time)) == 10


def test_bad_time():
    """Test raise exception if time stamp is unexpected."""
    badtime_filename = (
        _test_files_directory / "badtime_get_OBC_0_Data_1761872211189_1762025056586.csv"
    )
    with pytest.raises(ValueError):
        file_tools.read_file(badtime_filename)

    with pytest.raises(ValueError):
        file_tools.read_raw_file(badtime_filename)
