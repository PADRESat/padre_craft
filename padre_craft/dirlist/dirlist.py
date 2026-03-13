"""Functions to manage the directory or file list provided by the spacecraft"""

import re
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.table import QTable, Table
from astropy.time import Time
from astropy.timeseries import TimeSeries


class DirList:
    """
    Class to manage the directory or file list provided by the spacecraft.

    Parameters
    ----------
    dirlist_file: str or Path
        Path to the DirList file provided by the spacecraft. This file is usually an ASCII file containing a list of files currently stored on the on-board SD card, along with their sizes and timestamps.

    Example
    -------
    >>> from padre_craft.dirlist.dirlist import DirList
    >>> from padre_craft import _test_files_directory
    >>> dir_list = DirList(_test_files_directory / "padre_craft_dirlist_1772908542.txt")
    >>> print(len(dir_list))
    121
    >>> print(dir_list.file_count())
            name          count
            str23          int64
    ----------------------- -----
                    total   121
    padre_craft_padre_craft    27
            meddea_photon    32
                meddea_hk     4
            meddea_spectrum    33
                sharp_det0     1
                sharp_det1     1
                        ...   ...
                sharp_det5     1
                sharp_det6     9
                sharp_det7     1
            sharp_det_hk     2
            sharp_response     2
            sharp_histogram     1
        sharp_shipboot_hk     2
            sharp_ship_hk     1
    """

    def __init__(self, file_path: str | Path):
        if not isinstance(file_path, Path):
            file_path = Path(file_path)
        file_list = QTable(
            Table.read(
                file_path,
                format="ascii.csv",
                converters={
                    "size(in bytes)": int,
                    "file_name": str,
                    "timestamp": int,
                    "attributes": str,
                },
            )
        )
        file_list["size"] = (file_list["size(in bytes)"] * u.byte).to(u.MB)
        # filter out files with size = 0 bytes
        bool_array = file_list["size"] > 0 * u.MB
        file_list = file_list[bool_array]
        file_create_times = Time(file_list["timestamp"], format="unix", scale="utc")
        file_create_times.format = "isot"
        file_list["file_create_time"] = file_create_times
        file_list.meta["filename"] = str(file_path)

        match_short = re.search(r"(\d{10})", str(file_path.name))
        if match_short:
            time_unix_seconds = int(match_short.group(1))
            file_create_time = Time(time_unix_seconds, format="unix", scale="utc")
            file_list.meta["time"] = file_create_time.isot
        else:
            raise ValueError(
                f"Could not parse date from filename: {file_path}. "
                "Expected format: UNIX timestamp (10 digits)."
            )

        # normalize filenames - remove directory and "padre" and "_" from the filenames
        for i, this_f in enumerate(file_list["file_name"]):
            file_list["file_name"][i] = (
                this_f.replace("padre", "").replace("_", "").replace("/sd/", "")
            )

        self.file_list = file_list
        self.file_list["instrument"] = len(self.file_list) * ["padre_craft"]
        self.file_list["data_type"] = len(self.file_list) * ["padre_craft"]
        self.file_list["file_time"] = Time(
            ["2020-01-01T12:00:00Z"] * len(self.file_list), format="isot", scale="utc"
        )
        self._all_instr_data_types = {
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
        self._label_meddea_files()
        self._label_sharp_files()

    @classmethod
    def _parse_sharp_filename(cls, filename):
        if len(filename) == 21:
            m = re.match(
                r"^SP(\d{3}?)(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})\.(idx|dat)$",
                filename,
            )
        elif len(filename) == 20:
            m = re.match(
                r"^SP(\d{2}?)(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})\.(idx|dat)$",
                filename,
            )
        else:
            raise ValueError(f"Could not parse SHARP filename: {filename}")
        if m is not None:
            time_str = f"20{m[2]}-{m[3]}-{m[4]}T{m[5]}:{m[6]}:{m[7]}Z"
            APID = m[1]
            return APID, Time(time_str)
        else:
            raise ValueError(f"Could not parse SHARP filename: {filename}")

    @classmethod
    def _parse_meddea_filename(cls, filename):
        m = re.match(
            r"^MD(U8|A0|A2)(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})\.(dat)$",
            filename,
        )
        if m is not None:
            APID = filename[2:4]
            time_str = f"20{m[2]}-{m[3]}-{m[4]}T{m[5]}:{m[6]}:{m[7]}Z"
            return APID, Time(time_str)
        else:
            raise ValueError(f"Could not parse MeDDEA filename: {filename}")

    def __len__(self):
        return len(self.file_list)

    def _label_meddea_files(self):
        """Recognize and label all MeDDEA files based on the filename"""
        only_meddea_mask = np.array(
            [
                Path(this_f).name.startswith("MD")
                for this_f in self.file_list["file_name"]
            ]
        )
        only_idx_dat_mask = np.array(
            [
                Path(this_f).suffix in [".idx", ".dat"]
                for this_f in self.file_list["file_name"]
            ]
        )
        self.file_list["instrument"][only_meddea_mask & only_idx_dat_mask] = "meddea"
        for i, this_instrument in enumerate(self.file_list["instrument"]):
            filename = self.file_list["file_name"][i]
            if this_instrument == "meddea":
                this_apid, this_file_time = self._parse_meddea_filename(filename)
                self.file_list["file_time"][i] = this_file_time
                this_data_type = self._all_instr_data_types["meddea"][f"MD{this_apid}"]
                self.file_list["data_type"][i] = this_data_type

    def _label_sharp_files(self):
        """Recognize and label all SHARP files by filling in instrument and data type columns based on the filename"""
        only_sharp_mask = np.array(
            [
                Path(this_f).name.startswith("SP")
                for this_f in self.file_list["file_name"]
            ]
        )
        only_idx_dat_mask = np.array(
            [
                Path(this_f).suffix in [".idx", ".dat"]
                for this_f in self.file_list["file_name"]
            ]
        )
        self.file_list["instrument"][only_sharp_mask * only_idx_dat_mask] = "sharp"
        for i, this_instrument in enumerate(self.file_list["instrument"]):
            filename = self.file_list["file_name"][i]
            if this_instrument == "sharp":
                this_apid, this_file_time = self._parse_sharp_filename(filename)
                self.file_list["file_time"][i] = this_file_time
                this_data_type = self._all_instr_data_types["sharp"][f"SP{this_apid}"]
                self.file_list["data_type"][i] = this_data_type

    def available_instruments(self):
        return np.unique(self.file_list["instrument"])

    def available_data_types(self):
        return np.unique(self.file_list["data_type"])

    def _file_size_dict(self):
        result = {}
        result.update({"total": np.sum(self.file_list["size"])})
        for this_instrument, these_data_types in self._all_instr_data_types.items():
            for this_data_type in these_data_types.values():
                these_files = self.file_list[
                    (self.file_list["instrument"] == this_instrument)
                    & (self.file_list["data_type"] == this_data_type)
                ]
                total_size = np.sum(these_files["size"])
                result.update({f"{this_instrument}_{this_data_type}": total_size})
        return result

    def _file_count_dict(self):
        result = {}
        result.update({"total": len(self.file_list["size"])})
        for this_instrument, these_data_types in self._all_instr_data_types.items():
            for this_data_type in these_data_types.values():
                these_files = self.file_list[
                    (self.file_list["instrument"] == this_instrument)
                    & (self.file_list["data_type"] == this_data_type)
                ]
                result.update({f"{this_instrument}_{this_data_type}": len(these_files)})
        return result

    def file_size(self) -> QTable:
        file_size = self._file_size_dict()
        data = {"name": list(file_size.keys()), "size": list(file_size.values())}
        return QTable(data=data)

    def file_count(self) -> QTable:
        file_size = self._file_count_dict()
        data = {"name": list(file_size.keys()), "count": list(file_size.values())}
        return QTable(data=data)

    def __repr__(self):
        result = f"FileList {Path(self.file_list.meta['filename']).name} created on {self.file_list.meta['time']}.\n"
        result += f"Total size: {self._file_size_dict()['total']:.2f}\n"
        result += f"{self.file_size()}:\n"
        result += f"{self.file_count()}\n"
        result += f"{self.file_list}\n"
        return result

    def to_summary_ts(self, metric_type="size") -> TimeSeries:
        summary_ts = TimeSeries(time=[Time(self.file_list.meta["time"])])
        if metric_type == "size":
            data_dict = self._file_size_dict()
        elif metric_type == "count":
            data_dict = self._file_count_dict()
        else:
            raise ValueError(
                f"Invalid metric_type '{metric_type}'. Expected 'size' or 'count'."
            )
        for key, val in data_dict.items():
            if isinstance(val, u.Quantity):
                summary_ts[key] = val.value
            else:
                summary_ts[key] = val
        return summary_ts

    def only_sharp(self):
        """Return a new DirList object containing only SHARP files"""
        sharp_file_list = self.file_list[self.file_list["instrument"] == "sharp"]
        sharp_dirlist = DirList.__new__(DirList)
        sharp_dirlist.file_list = sharp_file_list
        sharp_dirlist._all_instr_data_types = self._all_instr_data_types
        return sharp_dirlist

    def only_meddea(self):
        """Return a new DirList object containing only MeDDEA files"""
        meddea_file_list = self.file_list[self.file_list["instrument"] == "meddea"]
        meddea_dirlist = DirList.__new__(DirList)
        meddea_dirlist.file_list = meddea_file_list
        meddea_dirlist._all_instr_data_types = self._all_instr_data_types
        return meddea_dirlist
