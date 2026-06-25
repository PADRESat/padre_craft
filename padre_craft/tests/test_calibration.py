from pathlib import Path
from unittest.mock import patch

import pytest
from astropy.io import fits
from astropy.table import Table

import padre_craft.calibration.calibration as calib
from padre_craft import _test_files_directory

test_file_paths = _test_files_directory.glob("padre_*.csv")


@pytest.mark.parametrize("this_path", list(test_file_paths))
def test_process_file(this_path, tmpdir, monkeypatch):
    # Set up the temporary directory as the current working directory
    monkeypatch.chdir(tmpdir)
    files = calib.process_file(this_path, output_fits=True)
    assert Path(files[0]).exists()
    # perform basic checks on the fits file
    hdul = fits.open(files[0])
    assert hdul[0].header["CREATOR"] == "padre_craft"
    # MEDDEA housekeeping data has extra rows
    if hdul[0].header["BTYPE"] == "meddea":
        assert len(hdul[1].data["timestamp_ms"]) == 13
    else:
        assert len(hdul[1].data["timestamp_ms"]) == 10


@patch("padre_craft.calibration.calibration.aws_db.record_dirlist")
def test_process_dirlist_file(mock_record_dirlist):
    """Test processing of dirlist files."""
    dirlist_file = "padre_craft_dirlist_1772908542.txt"
    test_dirlist_file = _test_files_directory / dirlist_file

    # Process the dirlist file
    output_files = calib.process_file(test_dirlist_file)

    # Check that the function returned a None placeholder
    assert len(output_files) == 1
    assert output_files[0] == "latest_dirlist.csv"

    file_list = Table.read(output_files[0], format="csv")
    assert Path(output_files[0]).exists()
    latest_file_list = Table.read(output_files[0], format="csv")
    assert len(latest_file_list) == 94
    assert set(latest_file_list['instrument']) == {'meddea', 'sharp'}


    # Verify that record_dirlist was called exactly once
    assert mock_record_dirlist.call_count == 1
