"""
A module for all things calibration.
"""

from pathlib import Path

from astropy.io import fits
from astropy.table import Table

import padre_craft.io.aws_db as aws_db
from padre_craft import log
from padre_craft.io import file_tools
from padre_craft.io.fits_tools import get_comment, get_obs_header, get_primary_header
from padre_craft.util.util import create_craft_filename

__all__ = [
    "process_file",
]


def process_file(filename: Path, overwrite=False) -> list:
    """
    This is the entry point for the pipeline processing.
    It runs all of the various processing steps required.

    Parameters
    ----------
    data_filename: str
        Fully specificied filename of an input file

    Returns
    -------
    output_filenames: list
        Fully specificied filenames for the output files.
    """
    log.info(f"Processing file {filename}.")

    output_files = []
    file_path = Path(filename)

    if file_path.suffix.lower() in [".csv"]:  # raw file
        data_ts = file_tools.read_file(file_path)

        # Prepare Metadata for Output Files Naming and Headers
        test_flag = False
        level_str = "l0"
        data_type = data_ts.meta["data_type"]

        aws_db.record_housekeeping(data_ts, data_type)
        data_table = Table(data_ts)

        # Get FITS Primary Header Template
        primary_hdr = get_primary_header(
            file_path, data_level=level_str, data_type=data_type
        )
        date_beg = data_ts.time[0]
        date_end = data_ts.time[-1]
        primary_hdr["DATE-BEG"] = (date_beg.fits, get_comment("DATE-BEG"))
        primary_hdr["DATE-END"] = (date_end.fits, get_comment("DATE-END"))
        primary_hdr["DATEREF"] = (date_beg.fits, get_comment("DATEREF"))

        colnames_to_remove = [
            "time",
        ]
        for this_col in colnames_to_remove:
            if this_col in data_table.colnames:
                data_table.remove_column(this_col)

        path = create_craft_filename(
            time=date_beg,
            level=level_str,
            descriptor=data_type,
            test=test_flag,
            version="1.0.0",
            overwrite=overwrite,
        )
        primary_hdr["FILENAME"] = (path, get_comment("FILENAME"))
        # record originating filename
        # aws_db.record_filename(file_path.name, date_beg, date_end)
        # record output filename
        # aws_db.record_filename(path.name, date_beg, date_end)
        empty_primary_hdu = fits.PrimaryHDU(header=primary_hdr)

        # Create HK HDU
        hk_header = get_obs_header(data_level=level_str, data_type=data_type)
        hk_header["DATE-BEG"] = (date_beg.fits, get_comment("DATE-BEG"))
        hk_header["DATEREF"] = (date_beg.fits, get_comment("DATEREF"))
        hk_header["FILENAME"] = (path, get_comment("FILENAME"))

        hk_hdu = fits.BinTableHDU(data=data_table, header=hk_header, name="HK")
        hk_hdu.add_checksum()
        hdul = fits.HDUList([empty_primary_hdu, hk_hdu])
        hdul.writeto(path, overwrite=overwrite, checksum=True)
        hdul.close()
        output_files.append(path)

    # add other tasks below
    return output_files
