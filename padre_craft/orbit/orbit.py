import csv
import urllib.request
from pathlib import Path

import astropy.units as u
import numpy as np
import skyfield.api
from astropy.time import Time
from astropy.timeseries import TimeSeries
from shapely import geometry
from skyfield.iokit import parse_tle_file

from padre_craft import NORAD_ID, _data_directory, log

__all__ = ["PadreOrbit", "get_ephemeris_file", "get_latest_tle", "get_celestrak_url"]

_SAA_VERTICES = np.array(
    [
        [-48.514111, -3.000365],
        [-79.502424, -16.25564],
        [-83.458183, -35.407033],
        [-69.116229, -43.917163],
        [-33.628958, -49.678745],
        [23.697636, -38.845567],
        [11.394957, -23.423524],
        [-48.514111, -3.000365],
    ]
)

_UPPER_BELT_VERTICES = np.array(
    [
        [-180, 66],
        [-120, 66],
        [-79, 62],
        [-20, 70],
        [32, 71],
        [112, 72],
        [180, 66],
        [180, 50],
        [150, 55],
        [85, 55],
        [45, 57],
        [-28, 56],
        [-100, 35],
        [-130, 41],
        [-180, 50],
        [-180, 66],
    ]
)

_LOWER_BELT_VERTICES = np.array(
    [
        [-180, -45],
        [-126, -48],
        [-83, -48],
        [-33, -49],
        [14, -48],
        [41, -51],
        [82, -45],
        [135, -40],
        [180, -45],
        [180, -60],
        [133, -54],
        [92, -61],
        [44, -72],
        [-29, -82],
        [-122, -82],
        [-180, -60],
        [-180, -45],
    ]
)


def vertices_to_polygon(vertices: np.array) -> geometry.Polygon:
    """Given a set of 2D vertices return a polygon.

    Parameters
    ----------
    vertices

    Returns
    -------
    polygon
    """
    # TODO: check that the polygon closes, the last point should be the same as first point
    # TODO: check that shape is 2
    pointList = []
    for this_lon, this_lat in vertices:
        this_point = geometry.Point(this_lon, this_lat)
        pointList.append(this_point)
    poly = geometry.Polygon(np.array(pointList))
    return poly


def get_ephemeris_file():
    """Get the path to the ephemeris file. If it does not exist, download it from JPL and put it into the data directory"""
    filename = "de421.bsp"
    file_path = _data_directory / filename
    if file_path.exists():
        return file_path
    else:
        url = f"https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/{filename}"
        log.info(f"{filename} is not found. Downloading from {url}")
        urllib.request.urlretrieve(url, filename)
        # move this file into the data directory
        new_path = _data_directory / filename
        Path(filename).replace(new_path)
        return new_path


def get_latest_tle() -> Path:
    """Get the latest TLE in CSV format. If it does not exist, download it from Celestrak and put it into the data directory.

    Returns
    -------
    tle_path
        Path to the tle file
    """
    url = get_celestrak_url(format="csv")
    filename = Path(f"{Time.now().iso.replace('-', '')[0:8]}_padre_tle.csv")
    file_path = _data_directory / filename
    if not file_path.exists():
        log.info(f"{filename} is not found. Downloading from {url}")
        urllib.request.urlretrieve(url, filename)
        Path(filename).replace(file_path)
    return file_path


def get_celestrak_url(format="csv") -> str:
    """Returns the url for the PADRE orbit information from celestrak

    Raises
    ------
    ValueError
        If format is not one recognized by Celestrak [JSON, CSV, TLE]

    Returns
    -------
    url: str
        The url to download the celestrak orbit information
    """
    if format.upper() not in ["CSV", "TLE", "JSON"]:
        raise ValueError(f"Format {format} is not recognized.")
    return f"https://celestrak.org/NORAD/elements/gp.php?CATNR={NORAD_ID}&FORMAT={format.upper()}"


class PadreOrbit:
    """A class to load and calculate the orbit of the PADRE spacecraft.

    Parameters
    ----------
    tle_filename : Path, str, or None
        Can parse CSV or standard TLE files. If given none, will download todays TLE from Celestrack

    Raises
    ------
    ValueError
        If does not recognize the TLE file or cannot download one if not given.

    Examples
    --------
    >>> from padre_craft.orbit import PadreOrbit
    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> from padre_craft import _test_files_directory
    >>> tle_filename = _test_files_directory / "20251219_padre_tle.csv"
    >>> padre_orbit = PadreOrbit(tle_filename)
    >>> padre_orbit.calculate(tstart=Time("2025-12-19T01:00"), tend=Time("2025-12-19T01:06"), dt=1 * u.min)
    >>> padre_orbit.in_sun
    array([ True,  True,  True,  True,  True,  True])
    >>> padre_orbit.in_particles
    array([False, False, False,  True,  True,  True])
    >>> padre_orbit.good_flag
    array([ True,  True,  True, False, False, False])
    """

    def __init__(self, tle_filename=None):
        self.ts = skyfield.api.load.timescale()
        self.satellite = None
        if tle_filename:
            if isinstance(tle_filename, str):
                tle_filename = Path(tle_filename)
        else:
            tle_filename = get_latest_tle()
        if tle_filename.suffix == ".csv":
            with skyfield.api.load.open(str(tle_filename), mode="r") as f:
                data = list(csv.DictReader(f))
            sats = [
                skyfield.api.EarthSatellite.from_omm(self.ts, fields) for fields in data
            ]
        elif tle_filename.suffix == ".tle":
            with skyfield.api.load.open(str(tle_filename)) as f:
                sats = list(parse_tle_file(f, self.ts))
        else:
            raise ValueError(f"File type of {tle_filename} not recognized.")
        if len(sats) >= 1:
            self.satellite = sats[0]
        else:
            raise ValueError(f"Problem loading tle file {tle_filename}")
        if self.satellite is None:
            raise ValueError("No TLE file found")

        self.eph = skyfield.api.load_file(get_ephemeris_file())
        self.polys = {
            "saa": vertices_to_polygon(_SAA_VERTICES),
            "upper_belt": vertices_to_polygon(_UPPER_BELT_VERTICES),
            "lower_belt": vertices_to_polygon(_LOWER_BELT_VERTICES),
        }

    def calculate(self, tstart: Time, tend: Time, dt=5 * u.s):
        """Calculates the position of PADRE in its orbit and associated parameters.

        Raises : ValueError
            If start or end times are beyond 2 weeks from TLE epoch
        """
        if (np.abs(tstart - self.satellite.epoch.to_astropy()) > 2 * u.week) or (
            np.abs(tend - self.satellite.epoch.to_astropy()) > 2 * u.week
        ):
            raise ValueError(
                f"Times provided too far from tle epoch {self.satellite.epoch}"
            )

        n_samples = np.ceil(((tend - tstart) / dt).decompose().value).astype(np.uint)
        self.timeseries = TimeSeries(
            time_start=tstart, time_delta=dt, n_samples=n_samples
        )
        t = self.ts.from_astropy(self.timeseries.time)
        self.times = self.timeseries.time
        self.geocentric = self.satellite.at(t)
        self.in_sun = self.satellite.at(t).is_sunlit(self.eph)
        geopos = skyfield.api.wgs84.geographic_position_of(self.geocentric)
        self.longitude = [
            float(this_lon) for this_lon in geopos.longitude.degrees
        ] * u.deg
        self.latitude = [
            float(this_lat) for this_lat in geopos.latitude.degrees
        ] * u.deg
        self.altitude = [float(this_lat) for this_lat in geopos.elevation.km] * u.km
        self.speed = (
            np.linalg.norm(self.geocentric.velocity.km_per_s, axis=0) * u.km / u.s
        )
        self.in_upper_belt = np.array(
            [
                geometry.Point(this_lon, this_lat).within(self.polys["upper_belt"])
                for this_lon, this_lat in zip(self.longitude.value, self.latitude.value)
            ]
        )
        self.in_lower_belt = np.array(
            [
                geometry.Point(this_lon, this_lat).within(self.polys["lower_belt"])
                for this_lon, this_lat in zip(self.longitude.value, self.latitude.value)
            ]
        )
        self.in_saa = np.array(
            [
                geometry.Point(this_lon, this_lat).within(self.polys["saa"])
                for this_lon, this_lat in zip(self.longitude.value, self.latitude.value)
            ]
        )
        self.in_particles = self.in_saa | self.in_lower_belt | self.in_upper_belt
        self.good_flag = self.in_sun & ~self.in_particles
        self.timeseries = TimeSeries(
            time=self.times,
            data={
                "longitude": self.longitude,
                "latitude": self.latitude,
                "altitude": self.altitude,
                "speed": self.speed,
                "in_sun": self.in_sun,
                "in_saa": self.in_saa,
                "in_upper_belt": self.in_upper_belt,
                "in_lower_belt": self.in_lower_belt,
            },
        )
