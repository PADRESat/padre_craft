import astropy.units as u
import pytest
from astropy.table import QTable
from astropy.time import Time

from padre_craft import _test_files_directory
from padre_craft.dirlist import dirlist

test_file = _test_files_directory / "padre_craft_dirlist_022426.txt"


@pytest.fixture
def file_list():
    return dirlist.read_dirlist(test_file)


def test_dirlist_class(file_list):
    assert isinstance(file_list, QTable)
    assert len(file_list) == 107
    assert file_list["size(in bytes)"].sum() > 840 * u.MB
    assert file_list.meta["filename"] == str(test_file)
    assert file_list.meta["time"] == "2026-02-24T00:00:00.000"


def test_meddealist_class(file_list):
    meddea_file_list = dirlist.MeDDEAFileList(file_list)
    assert isinstance(meddea_file_list, dirlist.MeDDEAFileList)
    assert len(meddea_file_list.file_list) == 69
    assert len(meddea_file_list.spec_files) == 33
    assert len(meddea_file_list.ph_files) == 32
    assert len(meddea_file_list.hk_files) == 4
    assert meddea_file_list.total_size() > 367 * u.MB
    assert meddea_file_list.total_size(data_type="spectrum") > 346 * u.MB
    assert meddea_file_list.total_size(data_type="photon") > 335 * u.MB
    assert meddea_file_list.total_size(data_type="hk") > 30 * u.MB


def test_parse_dirlist_filename_time():
    """Test the parse_dirlist_filename_time function with different formats."""
    # Test MMDDYY format (6 digits)
    time_6digit = dirlist.parse_dirlist_filename_time("Dirlist_Full_022426.txt")
    expected_6digit = Time("2026-02-24 00:00:00", format="iso", scale="utc")
    assert time_6digit == expected_6digit

    # Test MMDDYY format with different filename
    time_6digit2 = dirlist.parse_dirlist_filename_time("padre_craft_dirlist_022426.txt")
    assert time_6digit2 == expected_6digit

    # Test with Path object
    from pathlib import Path

    time_path = dirlist.parse_dirlist_filename_time(Path("test_dirlist_022426.txt"))
    assert time_path == expected_6digit

    # Test error case - no date pattern
    with pytest.raises(ValueError, match="Could not parse date from filename"):
        dirlist.parse_dirlist_filename_time("no_date_here.txt")


def test_summarize_dirlist(file_list):
    """Test the summarize_dirlist function."""
    summary = dirlist.summarize_dirlist(file_list)

    # Check that it returns a QTable
    assert isinstance(summary, QTable)

    # Check that it has the right columns
    assert "time" in summary.colnames
    assert "file_count_meddea" in summary.colnames
    assert "file_count_meddea_photon" in summary.colnames
    assert "file_count_meddea_spectrum" in summary.colnames
    assert "file_count_meddea_hk" in summary.colnames
    assert "file_count_sharp" in summary.colnames
    assert "file_count_total" in summary.colnames
    assert "size_meddea" in summary.colnames
    assert "size_meddea_photon" in summary.colnames
    assert "size_meddea_spectrum" in summary.colnames
    assert "size_meddea_hk" in summary.colnames
    assert "size_sharp" in summary.colnames
    assert "size_total" in summary.colnames

    # Check that it has one row
    assert len(summary) == 1

    # Check values
    assert summary["time"][0] == "2026-02-24T00:00:00.000"
    assert summary["file_count_total"][0] == 107
    assert summary["file_count_meddea"][0] == 69
    assert summary["file_count_meddea_photon"][0] > 0
    assert summary["file_count_meddea_spectrum"][0] > 0
    assert summary["file_count_meddea_hk"][0] > 0
    assert summary["file_count_sharp"][0] > 0  # There should be some SHARP files
    assert summary["size_total"][0] > 840
    assert summary["size_meddea"][0] > 0
    assert summary["size_meddea_photon"][0] > 0
    assert summary["size_meddea_spectrum"][0] > 0
    assert summary["size_meddea_hk"][0] > 0
    # Check that the sum of MeDDEA count components equals total MeDDEA count
    assert (
        summary["file_count_meddea_photon"][0]
        + summary["file_count_meddea_spectrum"][0]
        + summary["file_count_meddea_hk"][0]
    ) == summary["file_count_meddea"][0]
    # Check that the sum of MeDDEA components equals total MeDDEA size
    assert (
        summary["size_meddea_photon"][0]
        + summary["size_meddea_spectrum"][0]
        + summary["size_meddea_hk"][0]
    ) == summary["size_meddea"][0]
