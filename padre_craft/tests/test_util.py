import pytest

import padre_craft.util.util as util
from padre_craft import _test_files_directory
from padre_craft.io import read_file

TIME = "2024-04-06T12:06:21"
TIME_FORMATTED = "20240406T120621"

test_file_paths = _test_files_directory.glob("*.csv")
valid_datatypes = list(util.TOKEN_TO_DATATYPE.values())


def test_convert_meddea_colnames():
    test_file = (
        _test_files_directory
        / "padre_get_MEDDEA_HOUSE_KEEPING_Data_1762493454480_1762611866270.csv"
    )
    hk_ts = read_file(test_file)
    new_hk_ts = util.convert_meddea_colnames(hk_ts)
    assert "fp_temp" in new_hk_ts.colnames
    assert "dib_temp" in new_hk_ts.colnames
    assert "FPTemp" not in new_hk_ts.colnames
    assert "HVCurrent" not in new_hk_ts.colnames


@pytest.mark.parametrize("this_path", list(test_file_paths))
def test_filename_to_datatype(this_path):
    """Test that raw filenames output correct datatypes"""
    this_datatype = util.filename_to_datatype(this_path.name)
    assert isinstance(this_datatype, str)
    assert this_datatype in valid_datatypes


# fmt: off
@pytest.mark.parametrize("instrument,time,level,version,result", [
    ("meddea", TIME, "l1", "1.2.3", f"padre_meddea_l1_{TIME_FORMATTED}_v1.2.3.fits"),
    ("meddea", TIME, "l2", "2.4.5", f"padre_meddea_l2_{TIME_FORMATTED}_v2.4.5.fits"),
    ("sharp", TIME, "l2", "1.3.5", f"padre_sharp_l2_{TIME_FORMATTED}_v1.3.5.fits"),
    ("sharp", TIME, "l3", "2.4.5", f"padre_sharp_l3_{TIME_FORMATTED}_v2.4.5.fits"),
]
)
def test_science_filename_output_a(instrument, time, level, version, result):
    """Test simple cases with expected output.
    Since we are using the swxsoc create_science_filename, we are testing whether we did the config correctly in __init__.py"""
    assert (
        util.create_science_filename(instrument, time, level=level, version=version)
        == result
    )
# fmt: on


# fmt: off
@pytest.mark.parametrize("time,level,version,descriptor,result", [
    (TIME, "l1", "1.2.3", "adcs", f"padre_craft_l1_adcs_{TIME_FORMATTED}_v1.2.3.fits"),
    (TIME, "l2", "2.4.5", "adcs", f"padre_craft_l2_adcs_{TIME_FORMATTED}_v2.4.5.fits"),
    (TIME, "l2", "1.3.5", "gnss", f"padre_craft_l2_gnss_{TIME_FORMATTED}_v1.3.5.fits"),
    (TIME, "l3", "2.4.5", "hk", f"padre_craft_l3_hk_{TIME_FORMATTED}_v2.4.5.fits"),
]
)
def test_craft_filename_output_a(time, level, version, descriptor, result):
    """Test simple cases with expected output.
    Since we are using the swxsoc create_science_filename, we are testing whether we did the config correctly in __init__.py"""
    assert (
        util.create_craft_filename(time, level=level, descriptor=descriptor, version=version)
        == result
    )
# fmt: on
