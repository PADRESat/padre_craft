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
              time               longitude           latitude           altitude     ... in_sun in_saa in_upper_belt in_lower_belt
                                    deg                deg                 km        ...                                          
    ----------------------- ------------------- ------------------ ----------------- ... ------ ------ ------------- -------------
    2025-12-19T01:00:00.000 -152.34850804582567  72.93505067848139 603.3630125750256 ...   True  False         False         False
    2025-12-19T01:01:00.000 -157.47222559400785   69.5610001275033 602.8152817197313 ...   True  False         False         False
    2025-12-19T01:02:00.000 -161.25019407958612  66.08697774723791 602.1865865862586 ...   True  False         False         False
    2025-12-19T01:03:00.000  -164.1734992736378  62.55003870197339  601.489336473579 ...   True  False          True         False
    2025-12-19T01:04:00.000 -166.52663198202518  58.97083952464073 600.7375857366965 ...   True  False          True         False
    2025-12-19T01:05:00.000  -168.4832498821042 55.361670488292184  599.946794009871 ...   True  False          True         False

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
