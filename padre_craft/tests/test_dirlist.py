import astropy.units as u
import numpy as np
import pytest
from astropy.time import Time
from astropy.timeseries import TimeSeries

from padre_craft import _test_files_directory
from padre_craft.dirlist import dirlist

test_file = _test_files_directory / "padre_craft_dirlist_1772908542.txt"


@pytest.fixture
def dir_list():
    return dirlist.DirList(test_file)


def test_dirlist_class_input():
    dirlist.DirList(test_file)
    dirlist.DirList(str(test_file))


def test_dirlist_class(dir_list):
    assert isinstance(dir_list, dirlist.DirList)
    assert len(dir_list) == 107
    assert isinstance(dir_list.__repr__(), str)


def test_dirlist_contents(dir_list):
    assert dir_list._file_size_dict()["total"] > 840 * u.MB
    assert dir_list._file_size_dict()["total"] == np.sum(
        dir_list.file_list["size"]
    )
    assert dir_list.file_list.meta["filename"] == str(test_file)
    assert dir_list.file_list.meta["time"] == "2026-03-07T18:35:42.000"
    assert "file_create_time" in dir_list.file_list.colnames
    assert dir_list.file_list["file_create_time"][0] == Time(
        1769033135, format="unix", scale="utc"
    )  # Check that the first file's create time is correct


def test_dirlist_available_data_types(dir_list):
    data_types = dir_list.available_data_types()
    assert isinstance(data_types, np.ndarray)
    assert "photon" in data_types
    assert "spectrum" in data_types
    assert "hk" in data_types
    assert "160" in data_types
    assert "162" in data_types
    assert "padre_craft" in data_types


def test_dirlist_available_instruments(dir_list):
    instruments = dir_list.available_instruments()
    assert isinstance(instruments, np.ndarray)
    assert "meddea" in instruments
    assert "sharp" in instruments
    assert "padre_craft" in instruments


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
    [
        "total",
        "padre_craft_padre_craft",
        "meddea_photon",
        "meddea_hk",
        "meddea_spectrum",
        "sharp_162",
        "sharp_160",
    ],
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
    "name, count",
    [
        ("total", 107),
        ("padre_craft_padre_craft", 27),
        ("meddea_photon", 32),
        ("meddea_hk", 4),
        ("meddea_spectrum", 33),
        ("sharp_162", 10),
        ("sharp_160", 1),
    ],
)
def test_dirlist_file_count(dir_list, name, count):
    file_count_dict = dir_list._file_count_dict()
    #  file_count_table = dir_list.file_count()
    ts_count = dir_list.to_summary_ts(type="count")
    assert file_count_dict[name] == count
    assert ts_count[name] == count


@pytest.mark.parametrize(
    "name, size",
    [
        ("total", 861),
        ("padre_craft_padre_craft", 53),
        ("meddea_hk", 31),
        ("meddea_spectrum", 346),
        ("sharp_162", 94),
        ("sharp_160", 9.0e-05),
    ],
)
def test_dirlist_file_size(dir_list, name, size):
    file_size_dict = dir_list._file_size_dict()
    #  file_size_table = dir_list.file_size()
    ts_size = dir_list.to_summary_ts(type="size")
    assert file_size_dict[name] >= size * u.MB
    assert ts_size[name] >= size


def test_dirlist_to_summary_ts(dir_list):
    ts = dir_list.to_summary_ts(type="size")
    assert isinstance(ts, TimeSeries)
