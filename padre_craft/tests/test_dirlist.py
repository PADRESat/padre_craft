import astropy.units as u
import numpy as np
import pytest
from astropy.time import Time
from astropy.timeseries import TimeSeries

from padre_craft import _test_files_directory
from padre_craft.dirlist import dirlist

test_file = _test_files_directory / "padre_craft_dirlist_1772908542.txt"


_all_instr_data_types = {
    "padre_craft": {"padre_craft": "padre_craft"},
    "meddea": {"MDA0": "photon", "MDU8": "hk", "MDA2": "spectrum"},
    "sharp": {
        "SP10": "det0",
        "SP11": "det1",
        "SP12": "det2",
        "SP13": "det3",
        "SP14": "det4",
        "SP15": "det5",
        "SP16": "det6",
        "SP17": "det7",
        "SP20": "det_hk",
        "SP30": "response",
        "SP122": "histogram",
        "SP160": "shipboot_hk",
        "SP162": "ship_hk",
    },
}


@pytest.fixture
def dir_list():
    return dirlist.DirList(test_file)


def test_dirlist_class_input():
    dirlist.DirList(test_file)
    dirlist.DirList(str(test_file))


def test_dirlist_class(dir_list):
    assert isinstance(dir_list, dirlist.DirList)
    assert len(dir_list) == 121
    assert isinstance(dir_list.__repr__(), str)


def test_dirlist_contents(dir_list):
    assert dir_list._file_size_dict()["total"] > 840 * u.MB
    assert dir_list._file_size_dict()["total"] == np.sum(dir_list.file_list["size"])
    assert dir_list.file_list.meta["filename"] == str(test_file)
    assert dir_list.file_list.meta["time"] == "2026-03-07T18:35:42.000"
    assert "file_create_time" in dir_list.file_list.colnames
    assert dir_list.file_list["file_create_time"][0] == Time(
        1769033135, format="unix", scale="utc"
    )  # Check that the first file's create time is correct


def test_dirlist_available_data_types(dir_list):
    data_types = dir_list.available_data_types()
    assert isinstance(data_types, np.ndarray)
    for this_instrument in _all_instr_data_types.keys():
        for this_data_type in _all_instr_data_types[this_instrument].values():
            assert this_data_type in data_types


def test_dirlist_available_instruments(dir_list):
    instruments = dir_list.available_instruments()
    assert isinstance(instruments, np.ndarray)
    for this_instrument in _all_instr_data_types.keys():
        assert this_instrument in instruments


def test_dirlist_file_size_struct(dir_list):
    file_size_table = dir_list.file_size()
    assert isinstance(file_size_table, dirlist.QTable)
    assert "name" in file_size_table.colnames
    assert "size" in file_size_table.colnames
    assert len(file_size_table) > 0


def test_dirlist_file_count_struct(dir_list):
    file_count_table = dir_list.file_count()
    assert isinstance(file_count_table, dirlist.QTable)
    assert "name" in file_count_table.colnames
    assert "count" in file_count_table.colnames
    assert len(file_count_table) > 0


@pytest.mark.parametrize(
    "expected_row",
    ["total", "padre_craft_padre_craft"]
    + [f"sharp_{d}" for d in _all_instr_data_types["sharp"].values()]
    + [f"meddea_{d}" for d in _all_instr_data_types["meddea"].values()],
)
def test_dirlist_count_count_rows(dir_list, expected_row):
    file_count_dict = dir_list._file_count_dict()
    file_size_dict = dir_list._file_size_dict()
    file_count_table = dir_list.file_count()
    file_size_table = dir_list.file_size()
    assert expected_row in file_count_dict.keys()
    assert expected_row in file_size_dict.keys()
    assert expected_row in file_count_table["name"]
    assert expected_row in file_size_table["name"]


@pytest.mark.parametrize(
    "filename, expected",
    [
        (
            "SP160260205223638.idx",
            ("160", Time("2026-02-05T22:36:38.000", format="isot", scale="utc")),
        ),
        (
            "SP160270205223638.dat",
            ("160", Time("2027-02-05T22:36:38.000", format="isot", scale="utc")),
        ),
        (
            "SP16260205223639.dat",
            ("16", Time("2026-02-05T22:36:39.000", format="isot", scale="utc")),
        ),
        (
            "SP10260205224638.dat",
            ("10", Time("2026-02-05T22:46:38.000", format="isot", scale="utc")),
        ),
        (
            "SP160260205223938.dat",
            ("160", Time("2026-02-05T22:39:38.000", format="isot", scale="utc")),
        ),
        (
            "SP20260205223638.dat",
            ("20", Time("2026-02-05T22:36:38.000", format="isot", scale="utc")),
        ),
    ],
)
def test_parse_sharp_filename(filename, expected):
    result = dirlist.DirList._parse_sharp_filename(filename)
    assert result == expected


@pytest.mark.parametrize(
    "filename, expected",
    [
        (
            "MDA0260205223638.dat",
            ("A0", Time("2026-02-05T22:36:38.000", format="isot", scale="utc")),
        ),
        (
            "MDA2270205223638.dat",
            ("A2", Time("2027-02-05T22:36:38.000", format="isot", scale="utc")),
        ),
        (
            "MDU8260205223639.dat",
            ("U8", Time("2026-02-05T22:36:39.000", format="isot", scale="utc")),
        ),
        (
            "MDA0260205224638.dat",
            ("A0", Time("2026-02-05T22:46:38.000", format="isot", scale="utc")),
        ),
    ],
)
def test_parse_meddea_filename(filename, expected):
    result = dirlist.DirList._parse_meddea_filename(filename)
    assert result == expected


def test_parse_sharp_filename_error():
    filename = "test.csv"
    with pytest.raises(ValueError, match=f"Could not parse SHARP filename: {filename}"):
        dirlist.DirList._parse_sharp_filename(filename)


def test_parse_meddea_filename_error():
    filename = "test.csv"
    with pytest.raises(
        ValueError, match=f"Could not parse MeDDEA filename: {filename}"
    ):
        dirlist.DirList._parse_meddea_filename(filename)


@pytest.mark.parametrize(
    "name, count",
    [
        ("total", 121),
        ("padre_craft_padre_craft", 27),
        ("meddea_photon", 32),
        ("meddea_hk", 4),
        ("meddea_spectrum", 33),
        ("sharp_det0", 1),
        ("sharp_det1", 1),
        ("sharp_det2", 1),
        ("sharp_det3", 1),
        ("sharp_det4", 2),
        ("sharp_det5", 1),
        ("sharp_det6", 9),
        ("sharp_det7", 1),
        ("sharp_det_hk", 2),
        ("sharp_response", 2),
        ("sharp_histogram", 1),
        ("sharp_shipboot_hk", 2),
        ("sharp_ship_hk", 1),
    ],
)
def test_dirlist_file_count(dir_list, name, count):
    file_count_dict = dir_list._file_count_dict()
    #  file_count_table = dir_list.file_count()
    ts_count = dir_list.to_summary_ts(metric_type="count")
    assert file_count_dict[name] == count
    assert ts_count[name] == count


@pytest.mark.parametrize(
    "name, size",
    [
        ("total", 861),
        ("padre_craft_padre_craft", 53),
        ("meddea_hk", 31),
        ("meddea_spectrum", 346),
        ("sharp_det0", 10.49113),
        ("sharp_det1", 10.49113),
        ("sharp_det2", 10.49113),
        ("sharp_det3", 10.49113),
        ("sharp_det4", 14.41736),
        ("sharp_det5", 10.49113),
        ("sharp_det6", 94.388646),
        ("sharp_det7", 10.49113),
        ("sharp_det_hk", 12.087848),
        ("sharp_response", 10.491372),
        ("sharp_histogram", 10.49113),
        ("sharp_shipboot_hk", 10.491226),
        ("sharp_ship_hk", 3.2e-05),
    ],
)
def test_dirlist_file_size(dir_list, name, size):
    file_size_dict = dir_list._file_size_dict()
    #  file_size_table = dir_list.file_size()
    ts_size = dir_list.to_summary_ts(metric_type="size")
    assert file_size_dict[name] >= size * u.MB
    assert ts_size[name] >= size


def test_dirlist_to_summary_ts(dir_list):
    ts = dir_list.to_summary_ts(metric_type="size")
    assert isinstance(ts, TimeSeries)


def test_dirlist_to_summary_ts_invalid_metric(dir_list):
    with pytest.raises(
        ValueError, match="Invalid metric_type 'invalid'. Expected 'size' or 'count'."
    ):
        dir_list.to_summary_ts(metric_type="invalid")
