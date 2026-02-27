import astropy.units as u
import pytest
from astropy.table import QTable

from padre_craft import _test_files_directory
from padre_craft.dirlist import dirlist

test_file = _test_files_directory / "padre_craft_dirlist_022426.txt"


@pytest.fixture
def file_list():
    return dirlist.read_dirlist(test_file)


def test_dirlist_class(file_list):
    assert isinstance(file_list, QTable)
    assert len(file_list) == 106
    assert file_list["size(in bytes)"].sum() > 840 * u.MB


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
