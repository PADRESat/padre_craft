from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.time import Time
from astropy.timeseries import TimeSeries

import padre_craft.orbit.orbit as orbit
from padre_craft import _test_files_directory


@pytest.fixture
def padre_orbit():
    tle_filename = _test_files_directory / "20251219_padre_tle.csv"
    return orbit.PadreOrbit(tle_filename)


def test_orbit_class_withcsvfile(padre_orbit):
    assert isinstance(padre_orbit, orbit.PadreOrbit)
    assert padre_orbit.satellite is not None


def test_orbit_class_withtlefile(padre_orbit):
    tle_filename = _test_files_directory / "20251219_padre.tle"
    padre_orbit = orbit.PadreOrbit(tle_filename)
    assert isinstance(padre_orbit, orbit.PadreOrbit)
    assert padre_orbit.satellite is not None


def test_orbit_calculate_error(padre_orbit):
    # start time is bad
    with pytest.raises(ValueError):
        padre_orbit.calculate(
            tstart=Time("2025-11-19T00:00"), tend=Time("2025-12-19T01:00"), dt=1 * u.min
        )

    # end time is bad
    with pytest.raises(ValueError):
        padre_orbit.calculate(
            tstart=Time("2025-12-19T00:00"), tend=Time("2026-02-19T01:00"), dt=1 * u.min
        )


def test_orbit_calculate(padre_orbit):
    # 1 hour of data every 1 minute = 60 elements
    padre_orbit.calculate(
        tstart=Time("2025-12-19T00:00"), tend=Time("2025-12-19T01:00"), dt=1 * u.min
    )
    assert isinstance(padre_orbit.in_saa, np.ndarray)
    assert isinstance(padre_orbit.in_saa[0], np.bool)
    assert len(padre_orbit.in_saa) == 60
    assert len(padre_orbit.in_upper_belt) == 60
    assert len(padre_orbit.in_lower_belt) == 60
    assert len(padre_orbit.in_sun) == 60
    assert isinstance(padre_orbit.timeseries, TimeSeries)


def test_orbit_calculate2(padre_orbit):
    # 1 minute of data every 1 s = 60 elements
    padre_orbit.calculate(
        tstart=Time("2025-12-19T00:00"), tend=Time("2025-12-19T00:01"), dt=1 * u.s
    )
    assert isinstance(padre_orbit.in_saa, np.ndarray)
    assert isinstance(padre_orbit.in_saa[0], np.bool)
    assert len(padre_orbit.in_saa) == 60
    assert len(padre_orbit.in_upper_belt) == 60
    assert len(padre_orbit.in_lower_belt) == 60
    assert len(padre_orbit.in_sun) == 60
    assert isinstance(padre_orbit.timeseries, TimeSeries)


def test_orbit_calculate3(padre_orbit):
    """Test going across days"""
    # 1 day of data every 1 hour = 24 elements
    padre_orbit.calculate(
        tstart=Time("2025-12-19T00:00"), tend=Time("2025-12-20T00:00"), dt=1 * u.hr
    )
    assert isinstance(padre_orbit.in_saa, np.ndarray)
    assert isinstance(padre_orbit.in_saa[0], np.bool)
    assert len(padre_orbit.in_saa) == 24
    assert len(padre_orbit.in_upper_belt) == 24
    assert len(padre_orbit.in_lower_belt) == 24
    assert len(padre_orbit.in_sun) == 24
    assert isinstance(padre_orbit.timeseries, TimeSeries)


def test_orbit_calculate4(padre_orbit):
    """Check for physical quantities"""
    padre_orbit.calculate(
        tstart=Time("2025-12-19T00:00"), tend=Time("2025-12-19T00:05"), dt=1 * u.min
    )
    assert isinstance(padre_orbit.longitude, u.Quantity)
    assert padre_orbit.longitude.unit.is_equivalent(u.deg)
    assert isinstance(padre_orbit.latitude, u.Quantity)
    assert padre_orbit.latitude.unit.is_equivalent(u.deg)
    assert isinstance(padre_orbit.altitude, u.Quantity)
    assert padre_orbit.altitude.unit.is_equivalent(u.m)
    assert isinstance(padre_orbit.speed, u.Quantity)
    assert padre_orbit.speed.unit.is_equivalent(u.m / u.s)


@pytest.mark.remote_data
def test_orbit_class_nofile():
    padre_orbit = orbit.PadreOrbit()
    assert isinstance(padre_orbit, orbit.PadreOrbit)


@pytest.mark.remote_data
def test_get_latest_tle():
    tle_path = orbit.get_latest_tle()
    assert isinstance(tle_path, Path)
    assert tle_path.exists()


@pytest.mark.remote_data
def test_get_ephemeris_file():
    file_path = orbit.get_ephemeris_file()
    assert isinstance(file_path, Path)
    assert file_path.exists()


@pytest.mark.parametrize("format", ["json", "csv", "tle"])
def test_get_celestrak_url(format):
    url = orbit.get_celestrak_url(format=format)
    assert isinstance(url, str)
    assert url.count(str(orbit.NORAD_ID)) == 1


def test_get_celestrack_url_error():
    with pytest.raises(ValueError):
        orbit.get_celestrak_url(format="boo")
