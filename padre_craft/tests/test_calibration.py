from pathlib import Path

import pytest
from astropy.io import fits

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
