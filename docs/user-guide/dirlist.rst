.. _dirlist:

**************
DirList module
**************

Overview
========
The `~padre_craft.dirlist.dirlist` module provides tools to parse and analyze the directory listing files produced by the PADRE spacecraft.
These files are ASCII CSV files containing a list of files currently stored on the on-board SD card, along with their sizes and timestamps.
The module identifies files by instrument (MeDDEA, SHARP, or padre_craft) and data type, and exposes methods to summarize file counts and sizes.

The `~padre_craft.dirlist.dirlist.DirList` class is the primary interface for users to interact with the dirlist module.

.. code-block:: python

    >>> from padre_craft.dirlist import DirList

You can create an instance of the `~padre_craft.dirlist.dirlist.DirList` class by providing a dirlist file.
The filename must contain a 10-digit UNIX timestamp which is used to record when the listing was generated:

.. code-block:: python

    >>> from padre_craft import _test_files_directory
    >>> dirlist_file = _test_files_directory / "padre_craft_dirlist_1772908542.txt"
    >>> dir_list = DirList(dirlist_file)

The total number of files in the listing (excluding zero-byte files) can be obtained with :func:`len`:

.. code-block:: python

    >>> print(len(dir_list))
    121

The ``file_list`` attribute
===========================
The parsed file listing is stored in the ``file_list`` attribute as an `~astropy.table.QTable`.
Each row represents one file on the SD card. The table contains the following columns:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Column
     - Type
     - Description
   * - ``file_name``
     - str
     - Normalised filename (directory prefix, ``"padre"``, and underscores stripped).
   * - ``size(in bytes)``
     - int
     - Raw file size in bytes as reported by the spacecraft.
   * - ``size``
     - Quantity (MB)
     - File size converted to megabytes. Zero-byte files are excluded.
   * - ``timestamp``
     - int
     - File creation timestamp as a raw UNIX integer.
   * - ``file_create_time``
     - `~astropy.time.Time`
     - File creation time derived from ``timestamp``, in UTC.
   * - ``attributes``
     - str
     - Raw file attributes string from the dirlist.
   * - ``instrument``
     - str
     - Instrument that produced the file (``"meddea"``, ``"sharp"``, or ``"padre_craft"``).
   * - ``data_type``
     - str
     - Data product type within the instrument (e.g. ``"photon"``, ``"det0"``, ``"hk"``).
   * - ``file_time``
     - `~astropy.time.Time`
     - File creation time parsed from the filename itself (instrument-encoded timestamp).

You can access ``file_list`` directly to inspect or filter the full table:

.. code-block:: python

    >>> dir_list.file_list.colnames
    ['size(in bytes)', 'file_name', 'timestamp', 'attributes', 'size', 'file_create_time', 'instrument', 'data_type', 'file_time']

Because ``file_list`` is an `~astropy.table.QTable` you can use standard Astropy table indexing to
filter rows. For example, to list all MeDDEA photon files larger than 10 MB:

.. code-block:: python

    >>> mask = (dir_list.file_list["instrument"] == "meddea") & (dir_list.file_list["data_type"] == "photon") & (dir_list.file_list["size"] > 10 * u.MB)  # noqa: E501
    >>> print(dir_list.file_list[mask]["file_name", "size", "file_time"])
            file_name              size               file_time       
                            Mbyte
    -------------------- ------------------ -----------------------
    MDA0260121102925.dat 10.486877999999999 2026-01-21T10:29:25.000
    MDA0260121111943.dat 10.486673999999999 2026-01-21T11:19:43.000
    MDA0260121120821.dat          10.486794 2026-01-21T12:08:21.000
                    ...                ...                     ...
    MDA0260222192911.dat          10.486832 2026-02-22T19:29:11.000
    MDA0260223160137.dat 10.485885999999999 2026-02-23T16:01:37.000
    MDA0260209031209.dat          10.486296 2026-02-09T03:12:09.000
    MDA0260223171315.dat          10.486782 2026-02-23T17:13:15.000
    Length = 32 rows

Summarizing files
=================

File counts
-----------
The `~padre_craft.dirlist.dirlist.DirList.file_count` method returns a `~astropy.table.QTable` with the number of files
broken down by instrument and data type, as well as an overall total:

.. code-block:: python

    >>> dir_list.file_count().pprint(max_lines=-1)
    name          count
    ----------------------- -----
                    total   121
    padre_craft_padre_craft    27
            meddea_photon    32
                meddea_hk     4
            meddea_spectrum    33
                sharp_det0     1
                sharp_det1     1
                sharp_det2     1
                sharp_det3     1
                sharp_det4     2
                sharp_det5     1
                sharp_det6     9
                sharp_det7     1
            sharp_det_hk     2
            sharp_response     2
            sharp_histogram     1
        sharp_shipboot_hk     2
            sharp_ship_hk     1

File sizes
----------
The `~padre_craft.dirlist.dirlist.DirList.file_size` method returns a similar `~astropy.table.QTable` with the total
on-disk size (in MB) for each instrument and data type:

.. code-block:: python

    >>> dir_list.file_size().pprint(max_lines=-1)
                name                 size       
                                Mbyte
    ----------------------- ------------------
                    total         982.645443
    padre_craft_padre_craft  53.62375699999999
            meddea_photon 335.56548599999996
                meddea_hk 31.637251999999997
            meddea_spectrum         346.504554
                sharp_det0           10.49113
                sharp_det1           10.49113
                sharp_det2           10.49113
                sharp_det3           10.49113
                sharp_det4           14.41736
                sharp_det5           10.49113
                sharp_det6          94.388646
                sharp_det7           10.49113
            sharp_det_hk          12.087848
            sharp_response          10.491372
            sharp_histogram           10.49113
        sharp_shipboot_hk          10.491226
            sharp_ship_hk            3.2e-05

Available instruments and data types
-------------------------------------
You can query which instruments and data types are present in a given listing:

.. code-block:: python

    >>> print(dir_list.available_instruments())
     instrument
    -----------
        meddea
    padre_craft
        sharp
    >>> print(dir_list.available_data_types())
    ['det0' 'det1' 'det2' 'det3' 'det4' 'det5' 'det6' 'det7' 'det_hk'
     'histogram' 'hk' 'padre_craft' 'photon' 'response' 'ship_hk'
     'shipboot_hk' 'spectrum']

Filtering by instrument
=======================
The `~padre_craft.dirlist.dirlist.DirList.only_sharp` and `~padre_craft.dirlist.dirlist.DirList.only_meddea` methods
return a new `~padre_craft.dirlist.dirlist.DirList` instance containing only files from the specified instrument:

.. code-block:: python

    >>> sharp_files = dir_list.only_sharp()
    >>> print(len(sharp_files))
    26
    >>> meddea_files = dir_list.only_meddea()
    >>> print(len(meddea_files))
    69

Converting to a TimeSeries
==========================
The `~padre_craft.dirlist.dirlist.DirList.to_summary_ts` method converts the summary statistics to an
`~astropy.timeseries.TimeSeries` object, which can be uploaded to a time-series database.
The ``metric_type`` parameter selects between ``"size"`` (default) and ``"count"``:

.. code-block:: python

    >>> ts = dir_list.to_summary_ts(metric_type="count")
    >>> print(ts)
    <TimeSeries length=1>
              time             total padre_craft_padre_craft ... sharp_ship_hk
    ----------------------- ------- ----------------------- ... -------------
    2026-03-13T01:22:22.000     121                      27 ...             1
