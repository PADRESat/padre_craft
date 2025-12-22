.. _orbit:

************
Orbit module
************

Overview
========
The `~padre_craft.orbit.orbit` module provides tools to compute and analyze the orbit of the Padre Craft satellite.
It includes functionalities to calculate the satellite's position, velocity, and various orbital parameters over time.
The module leverages `skyfield` for precise orbital mechanics calculations.
It bases its orbit determination on Two-Line Element (TLE) data obtained from Celestrak.
To determine whether the satellite is in sunlight, it makes use of ephermis data provided by JPL.

The `~padre_craft.orbit.orbit.PadreOrbit` class is the primary interface for users to interact with the orbit module.

.. code-block:: python

    >>> from padre_craft.orbit import PadreOrbit

You can create an instance of the `~padre_craft.orbit.orbit.PadreOrbit` class by providing a TLE file:

.. code-block:: python

    >>> from padre_craft import _test_files_directory
    >>> tle_filename = _test_files_directory / "20251219_padre_tle.csv"
    >>> padre_orbit = PadreOrbit(tle_filename)

Once you have an instance of `PadreOrbit`, you can calculate the orbit over a specified time range:

.. code-block:: python

    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> padre_orbit.calculate(tstart=Time("2025-12-19T01:00"), tend=Time("2025-12-19T01:06"), dt=1 * u.min)

Once the orbit is calculated, you can access various properties of the orbit, such as position, velocity, and flags indicating whether the satellite is in sunlight or within radiation belts.

.. code-block:: python

    >>> padre_orbit.in_sun
    array([ True,  True,  True,  True,  True,  True])
    >>> padre_orbit.in_particles
    array([False, False, False,  True,  True,  True])
    >>> padre_orbit.good_flag
    array([ True,  True,  True, False, False, False])

Most state information is stored as boolean arrays, where each element corresponds to a time step in the calculated orbit.
A timeseries of the orbit data can be obtained using the `timeseries` parameter:

.. code-block:: python

    >>> padre_orbit_ts = padre_orbit.timeseries
    >>> print(padre_orbit_ts)
              time               longitude      ... in_upper_belt in_lower_belt
                                    deg         ...
    ----------------------- ------------------- ... ------------- -------------
    2025-12-19T01:00:00.000 -152.34855524057318 ...         False         False
    2025-12-19T01:01:00.000 -157.47227278912428 ...         False         False
    2025-12-19T01:02:00.000  -161.2502412750741 ...         False         False
    2025-12-19T01:03:00.000 -164.17354646949326 ...          True         False
    2025-12-19T01:04:00.000 -166.52667917825477 ...          True         False
    2025-12-19T01:05:00.000 -168.48329707870266 ...          True         False

Visualization
=============
The `PadreOrbit` class includes a method to visualize the orbit state over time.
We will calculate a larger orbit time range in the following examples.

You can plot the orbit state using the `plot_state` method:

.. plot::
    :include-source:

    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> from padre_craft.orbit import PadreOrbit
    >>> from padre_craft import _test_files_directory
    >>> tle_filename = _test_files_directory / "20251219_padre_tle.csv"
    >>> padre_orbit = PadreOrbit(tle_filename)
    >>> padre_orbit.calculate(tstart=Time("2025-12-19T02:45"), tend=Time("2025-12-19T03:45"), dt=1 * u.min)
    >>> padre_orbit.plot_state()

This will generate a series of subplots showing the satellite's status regarding sunlight exposure, radiation belt crossings, and overall data quality over the specified time range.

You can also visualize the satellite's geolocation on a world map using the `plot_geolocation` method:

.. plot::
    :include-source:

    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> from padre_craft.orbit import PadreOrbit
    >>> from padre_craft import _test_files_directory
    >>> tle_filename = _test_files_directory / "20251219_padre_tle.csv"
    >>> padre_orbit = PadreOrbit(tle_filename)
    >>> padre_orbit.calculate(tstart=Time("2025-12-19T02:45"), tend=Time("2025-12-19T03:45"), dt=1 * u.min)
    >>> padre_orbit.plot_geolocation()
