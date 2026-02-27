"""Functions to manage the directory or file list provided by the spacecraft"""

import re
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.table import QTable, Table, vstack
from astropy.time import Time, TimeDelta
from astropy.timeseries import TimeSeries

from padre_craft.orbit import PadreOrbit


def read_dirlist(dirlist_file):
    """
    Reads a DirList file and returns a QTable with the file information.

    Parameters
    ---
    dirlist_file: str or Path
        Path to the DirList file provided by the spacecraft. This file is usually an ASCII file containing a list of files currently stored on the on-board SD card, along with their sizes and timestamps.

    Returns
    ---
    file_list: QTable
        Table containing the file information from the DirList, including "size(in bytes)", "file_name", "timestamp", and "attributes" columns.
    """
    file_list = QTable(
        Table.read(
            dirlist_file,
            format="ascii.csv",
            converters={
                "size(in bytes)": int,
                "file_name": str,
                "timestamp": int,
                "attributes": str,
            },
        )
    )
    file_list["size(in bytes)"] = (file_list["size(in bytes)"] * u.byte).to(u.MB)
    # filter out files with size = 0 bytes
    bool_array = file_list["size(in bytes)"] > 0 * u.MB
    file_list = file_list[bool_array]
    return file_list


class MeDDEAFileList:
    """
    Class to manage the directory or file list provided by the spacecraft.
    This class provides methods to read in the DirList, parse the filenames to extract timestamps, and add orbit information based on the timestamps.

    Parameters
    ----------
    dirlist_file: str or Path
        Path to the DirList file provided by the spacecraft. This file is usually an ASCII file containing a list of files currently stored on the on-board SD card, along with their sizes and timestamps.

    Example
    -------

    """

    def __init__(self, file_list: QTable):
        # Keep only the MeDDEA files:
        only_meddea_mask = np.array(
            ["MD" in Path(this_f).name for this_f in file_list["file_name"]]
        )
        no_idx_mask = np.array(
            ["dat" in Path(this_f).name for this_f in file_list["file_name"]]
        )
        meddea_file_list = file_list[only_meddea_mask * no_idx_mask]
        # remove the directory from the filenames
        meddea_file_list["file_name"] = [
            Path(this_f).name for this_f in meddea_file_list["file_name"]
        ]
        # normalize filenames
        for i, this_f in enumerate(meddea_file_list["file_name"]):
            meddea_file_list["file_name"][i] = this_f.replace("padre", "").replace(
                "_", ""
            )

        file_type = []
        for this_row in meddea_file_list:
            if "U8" in this_row["file_name"]:
                file_type.append("hk")
            elif "A0" in this_row["file_name"]:
                file_type.append("photon")
            elif "A2" in this_row["file_name"]:
                file_type.append("spectrum")

        meddea_file_list["data_type"] = file_type
        start_time = [
            parse_science_filename(Path(this_row["file_name"]).name)["time"]
            for this_row in meddea_file_list
        ]
        meddea_files = TimeSeries(data=meddea_file_list, time=start_time)
        meddea_files.sort()

        self.file_list = meddea_files

    def __repr__(self):
        result = f"MeDDEAFileList with {len(self.file_list)} files, including {len(self.hk_files)} housekeeping files, {len(self.ph_files)} photon files, and {len(self.spec_files)} spectrum files.\n"
        tdiff = self.file_list["time"][-1] - self.file_list["time"][0]
        result += f"Time range: {self.file_list['time'][0].iso} to {self.file_list['time'][-1].iso} ({tdiff.to(u.day):.2f} days)\n"
        tdiff = self.file_list["time"][-1] - self.file_list["time"][0]
        result += "Size for each data type:\n"
        result += f"  Housekeeping: {self.total_size('hk'):.2f}\n"
        result += f"  Photon: {self.total_size('photon'):.2f}\n"
        result += f"  Spectrum: {self.total_size('spectrum'):.2f}\n"
        result += f"Total size: {self.total_size():.2f}\n"
        result += f"Files:\n{self.file_list}"
        return result

    @property
    def hk_files(self):
        return self.file_list[self.file_list["data_type"] == "hk"]

    @property
    def ph_files(self):
        return self.file_list[self.file_list["data_type"] == "photon"]

    @property
    def spec_files(self):
        return self.file_list[self.file_list["data_type"] == "spectrum"]

    def total_size(self, data_type=None):
        if data_type is not None:
            return np.sum(
                self.file_list[self.file_list["data_type"] == data_type][
                    "size(in bytes)"
                ]
            )
        else:
            return np.sum(self.file_list["size(in bytes)"])

    @classmethod
    def parse_file_name(cls, filename):
        """
        Parse MeDDEA science filenames from the DirList and extract their timestamps.
        The DirList often contains both "padreMDA*_YYMMDDHHMMSS.dat" and "MDA*YYMMDDHHMMSS.dat" formats.

        Parameters
        ---
        filename: str
            Filename to parse.

        Returns
        ---
        time_dict: dict
            Dictionary with "time" key containing astropy Time object.
        """

        pattern1 = re.search(
            r"padreMD[UA][0-9]_(\d{12})\.dat", filename
        )  # encompasses A0, A2, U8
        pattern2 = re.search(
            r"MD[UA][0-9](\d{12})\.dat", filename
        )  # encompasses A0, A2, U8

        if pattern1:
            timestamp_str = pattern1.group(1)
        elif pattern2:
            timestamp_str = pattern2.group(1)
        else:
            raise ValueError(f"Filename {filename} does not match expected patterns")

        # Parse timestamp: YYMMDDHHMMSS
        year = "20" + timestamp_str[0:2]
        month = timestamp_str[2:4]
        day = timestamp_str[4:6]
        hour = timestamp_str[6:8]
        minute = timestamp_str[8:10]
        second = timestamp_str[10:12]

        time_str = f"{year}-{month}-{day} {hour}:{minute}:{second}"
        time_obj = Time(time_str, format="iso", scale="utc")

        time_dict = {"time": time_obj}

        return time_dict

    def add_durations(self, remove_suspect_durations=True):
        """
        Add end_time and duration columns to a list of files with time information.
        Remove the last file from the list since it does not have an end time or duration.
        """
        data_types = np.unique(self.file_list["data_type"])
        my_lists = []
        for this_data_type in data_types:
            this_list = self.file_list[
                self.file_list["data_type"] == this_data_type
            ].copy()

            if len(this_list) == 1:
                print(
                    "List must contain more than one file to calculate end times and durations."
                )
            else:
                end_times = self.file_list["time"][1:]
                durations = (
                    self.file_list["time"][1:] - self.file_list["time"][:-1]
                ).to("min")
                new_list = self.file_list[:-1]

                new_list["end_time"] = end_times
                new_list["duration"] = durations
                my_lists.append(new_list)

        self.file_list = vstack(my_lists)
        print(self.file_list)
        if remove_suspect_durations:
            valid_durations = {"spectrum": 71 * u.min, "photon": 150 * u.min}
            bool_array = np.ones(len(self.file_list), dtype=bool)
            for data_type, max_duration in valid_durations.items():
                my_files_bool = self.file_list[self.file_list["data_type"] == data_type]
                bad_duration_bool = self.file_list["duration"] > max_duration
                if any(bad_duration_bool * my_files_bool):
                    bad_files_bool = bad_duration_bool * my_files_bool
                    print(
                        f"Warning, {data_type} should be {max_duration} long and not longer"
                    )
                    bool_array *= ~bad_files_bool
            self.file_list = self.file_list[bool_array]


class SHARPFileList:
    """
    Class to manage the directory or file list provided by the spacecraft, specifically for SHARP files.
    This class provides methods to read in the DirList, parse the filenames to extract timestamps, and add orbit information based on the timestamps.

    Parameters
    ----------
    file_list: QTable
        Path to the DirList file provided by the spacecraft. This file is usually an ASCII file containing a list of files currently stored on the on-board SD card, along with their sizes and timestamps.

    Example
    -------

    """

    def __init__(self, file_list: QTable):
        # Keep only the SHARP files:
        only_sharp_mask = np.array(
            ["SP1" in Path(this_f).name for this_f in file_list["file_name"]]
        )
        sharp_file_list = file_list[only_sharp_mask]
        # remove the directory from the filenames
        sharp_file_list["file_name"] = [
            Path(this_f).name for this_f in sharp_file_list["file_name"]
        ]
        self.file_list = sharp_file_list

    def total_size(self):
        return np.sum(self.file_list["size(in bytes)"])


def parse_science_filename(filename):
    """
    Parse MeDDEA science filenames from the DirList and extract their timestamps.
    The DirList often contains both "padreMDA*_YYMMDDHHMMSS.dat" and "MDA*YYMMDDHHMMSS.dat" formats.

    Parameters
    ---
    filename: str
        Filename to parse.

    Returns
    ---
    time_dict: dict
        Dictionary with "time" key containing astropy Time object.
    """

    pattern1 = re.search(
        r"padreMD[UA][0-9]_(\d{12})\.dat", filename
    )  # encompasses A0, A2, U8
    pattern2 = re.search(
        r"MD[UA][0-9](\d{12})\.dat", filename
    )  # encompasses A0, A2, U8

    if pattern1:
        timestamp_str = pattern1.group(1)
    elif pattern2:
        timestamp_str = pattern2.group(1)
    else:
        raise ValueError(f"Filename {filename} does not match expected patterns")

    # Parse timestamp: YYMMDDHHMMSS
    year = "20" + timestamp_str[0:2]
    month = timestamp_str[2:4]
    day = timestamp_str[4:6]
    hour = timestamp_str[6:8]
    minute = timestamp_str[8:10]
    second = timestamp_str[10:12]

    time_str = f"{year}-{month}-{day} {hour}:{minute}:{second}"
    time_obj = Time(time_str, format="iso", scale="utc")

    time_dict = {"time": time_obj}

    return time_dict


def add_end_durations(this_list):
    """
    Add end_time and duration columns to a list of files with time information.
    Remove the last file from the list since it does not have an end time or duration.

    Parameters
    ---
    this_list: QTable
        Table containing at least a "time" column with astropy Time objects.

    Returns
    ---
    this_list: QTable
        Updated table with "end_time" and "duration" columns added.
    """
    if len(this_list) > 1:
        new_list = this_list.copy()
        end_times = this_list["time"][1:]
        durations = (this_list["time"][1:] - this_list["time"][:-1]).to("min")

        new_list = this_list[:-1]

        new_list["end_time"] = end_times
        new_list["duration"] = durations
        return new_list
    else:
        print(
            "List must contain more than one file to calculate end times and durations."
        )
        return this_list


def parse_dirlist(dirlist_file, old_filter=2 * u.week, remove_suspect_durations=True):
    """
    Reads a DirList and returns sorted lists of MeDDEA housekeeping, photon, and spectrum files.

    Parameters
    ---
    dirlist: ascii file
        ASCII file provided by the PADRE operator (usually Ace), containing a list of files currently stored on the on-board SD card.

    Returns
    ---
    hk_files:
        Sorted list of housekeeping files contained in the DirList.

    ph_files:
        Sorted list of photon files contained in the DirList.

    spec_files:
        Sorted list of spectrum files contained in the DirList.
    """
    # Read in the DirList:
    file_list = QTable(
        Table.read(
            dirlist_file,
            format="ascii.csv",
            converters={
                "size(in bytes)": int,
                "file_name": str,
                "timestamp": int,
                "attributes": str,
            },
        )
    )
    file_list["size(in bytes)"] = (file_list["size(in bytes)"] * u.byte).to(u.MB)
    # filter out files with size = 0 bytes
    bool_array = file_list["size(in bytes)"] > 0 * u.MB
    file_list = file_list[bool_array]

    # Keep only the MeDDEA files:
    only_meddea_mask = np.array(
        ["MD" in Path(this_f).name for this_f in file_list["file_name"]]
    )
    no_idx_mask = np.array(
        ["dat" in Path(this_f).name for this_f in file_list["file_name"]]
    )
    meddea_file_list = file_list[only_meddea_mask * no_idx_mask]
    # remove the directory from the filenames
    meddea_file_list["file_name"] = [
        Path(this_f).name for this_f in meddea_file_list["file_name"]
    ]
    # normalize filenames
    for i, this_f in enumerate(meddea_file_list["file_name"]):
        meddea_file_list["file_name"][i] = this_f.replace("padre", "").replace("_", "")

    file_type = []
    for this_row in meddea_file_list:
        if "U8" in this_row["file_name"]:
            file_type.append("hk")
        elif "A0" in this_row["file_name"]:
            file_type.append("photon")
        elif "A2" in this_row["file_name"]:
            file_type.append("spectrum")

    meddea_file_list["data_type"] = file_type
    start_time = [
        parse_science_filename(Path(this_row["file_name"]).name)["time"]
        for this_row in meddea_file_list
    ]
    meddea_files = TimeSeries(data=meddea_file_list, time=start_time)
    meddea_files.sort()
    if old_filter is not None:
        old_time = Time.now() - old_filter
        meddea_files = meddea_files[meddea_files["time"] > old_time]

    hk_files = add_end_durations(meddea_files[meddea_files["data_type"] == "hk"])
    ph_files = add_end_durations(meddea_files[meddea_files["data_type"] == "photon"])
    spec_files = add_end_durations(
        meddea_files[meddea_files["data_type"] == "spectrum"]
    )

    if any(spec_files["duration"] > 71 * u.min):
        bad_durations = spec_files["duration"] > 72 * u.min
        print(spec_files[bad_durations])
        print("Warning, spec should should be 71 minutes long and not longer")
        if remove_suspect_durations:
            spec_files = spec_files[~bad_durations]
    if any(ph_files["duration"] > 150 * u.min):
        bad_durations = ph_files["duration"] > 150 * u.min
        print(ph_files[bad_durations])
        print("Warning, ph should should be 150 minutes long and not longer")
        if remove_suspect_durations:
            ph_files = ph_files[~bad_durations]

    return hk_files, ph_files, spec_files


def add_orbit_info(file_list):
    """
    Add orbit information to a list of MeDDEA files based on their timestamps.
    This function uses the PadreOrbit class to determine the orbit number for each file based on its timestamp.

    Parameters
    ---
    file_list: QTable
        Table containing at least a "time" column with astropy Time objects, and a "file_name" column with the corresponding filenames.

    Returns
    ---
    file_list: QTable
        Updated table with percent of good sun observations and good calibration observations added as "good_sun_obs" and "good_cal_obs" columns, respectively.
    """
    # evaluate photons files for calibration files
    new_file_list = file_list.copy()
    padre_orbit = PadreOrbit()
    good_sun_obs_list = [0.0] * len(file_list)
    good_cal_obs_list = [0.0] * len(file_list)
    for i, this_row in enumerate(file_list):
        if this_row["time"] >= (Time.now() - TimeDelta(10 * u.day)):
            padre_orbit.calculate(tstart=this_row["time"], tend=this_row["end_time"])
            ts = padre_orbit.timeseries
            in_particles = ts["in_saa"] | ts["in_upper_belt"] | ts["in_lower_belt"]
            good_sun_obs = ts["in_sun"] * (~in_particles)
            good_cal_obs = ~ts["in_sun"] * (~in_particles)
            good_sun_obs_list[i] = np.sum(good_sun_obs) / len(ts)
            good_cal_obs_list[i] = np.sum(good_cal_obs) / len(ts)

            # calculate the percent of good data in the photon file based on the percentage of time spent in eclipse vs. in sunlight, since we only want to keep photon files that are mostly in sunlight for calibration purposes
    new_file_list["good_sun_obs"] = good_sun_obs_list * u.percent * 100
    new_file_list["good_cal_obs"] = good_cal_obs_list * u.percent * 100
    return new_file_list
