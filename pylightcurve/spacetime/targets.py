
__all__ = ['FixedTarget']

import numpy as np

from pylightcurve.errors import *
from pylightcurve.__databases__ import plc_data
from pylightcurve.analysis.curve_fit import curve_fit
from pylightcurve.spacetime.angles import arccos, Rad, _request_angle


class _Target:

    def __init__(self, ra, dec):

        _request_angle(ra)
        _request_angle(dec)

        if 270 > dec.deg() > 90:
            raise PyLCInputError('Declination must within -90, 90 degrees')

        self.dec = dec
        self.ra = ra

        self.coord = '{0} {1}'.format(self.ra.hms(), self.dec.dms_coord())


class FixedTarget(_Target):
    """A fixed position on the sky, for time conversions.

    The time conversions need the position of the target, because the light-travel
    time to the barycentre or to the centre of the Sun depends on the direction.

    Parameters
    ----------
    ra : Degrees or Hours or Rad
        Right ascension.
    dec : Degrees or Hours or Rad
        Declination, between -90 and 90 degrees.

    Attributes
    ----------
    ra, dec : Degrees or Hours or Rad
        The position.
    coord : str
        The position as ``'hh:mm:ss.s +dd:mm:ss.s'``.

    Raises
    ------
    PyLCInputError
        If ``ra`` or ``dec`` is not an angle object, or if ``dec`` is outside -90 to 90
        degrees.

    Notes
    -----
    The time systems are:

    - ``JD_UTC``, ``MJD_UTC``: Julian Date or Modified Julian Date (JD - 2400000.5) at
      the observer, in UTC.
    - ``HJD_UTC``, ``HJD_TDB``: Heliocentric Julian Date, the time at the centre of the
      Sun, in UTC or TDB.
    - ``BJD_UTC``, ``BJD_TDB``: Barycentric Julian Date, the time at the barycentre of
      the Solar System, in UTC or TDB. BJD_TDB is the system used internally by
      PyLightcurve.
    - ``HJD_TT`` and ``BJD_TT`` are treated as ``HJD_TDB`` and ``BJD_TDB``.

    The observer is assumed to be at the centre of the Earth.

    Examples
    --------
    >>> import pylightcurve as plc
    >>> target = plc.FixedTarget(plc.Hours('22:03:10.7729'), plc.Degrees('+18:53:03.548'))
    >>> target.convert_to_bjd_tdb(2458485.00380255, 'HJD_UTC')

    """

    def __init__(self, ra, dec):

        _Target.__init__(self, ra, dec)

    def convert_to_bjd_tdb(self, time, time_format):
        """Convert times to BJD_TDB.

        Parameters
        ----------
        time : float or array_like
            Time(s) to convert.
        time_format : str
            Time system of ``time``: ``'BJD_TDB'``, ``'BJD_TT'``,
            ``'BJD_UTC'``, ``'HJD_TDB'``, ``'HJD_TT'``, ``'HJD_UTC'``, ``'JD_UTC'`` or
            ``'MJD_UTC'``.

        Returns
        -------
        float or numpy.ndarray
            The time(s) in BJD_TDB: a float for a float input, an array otherwise.

        Raises
        ------
        PyLCInputError
            If ``time_format`` is not valid or ``time`` cannot be converted to floats.

        """

        if isinstance(time, float):

            if time_format in ['BJD_TDB', 'BJD_TT']:
                return time
            elif time_format == 'JD_UTC':
                return self._bjd_tdb(None, time)
            elif time_format == 'MJD_UTC':
                return self._bjd_tdb(None, time + 2400000.5)
            elif time_format == 'BJD_UTC':
                return self._bjd_tdb(None, curve_fit(self._bjd_utc, [0], [time], p0=[time])[0][0])
            elif time_format in ['HJD_TDB', 'HJD_TT']:
                return self._bjd_tdb(None, curve_fit(self._hjd_tdb, [0], [time], p0=[time])[0][0])
            elif time_format == 'HJD_UTC':
                return self._bjd_tdb(None, curve_fit(self._hjd_utc, [0], [time], p0=[time])[0][0])
            else:
                raise PyLCInputError(
                    'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                    'HJD_UTC, HJD_BJD, HJD_TT')

        else:
            try:
                time = np.array(time, dtype=float)

                if time_format in ['BJD_TDB', 'BJD_TT']:
                    return time
                elif time_format == 'JD_UTC':
                    return np.array([self._bjd_tdb(None, ff) for ff in time])
                elif time_format == 'MJD_UTC':
                    return np.array([self._bjd_tdb(None, ff + 2400000.5) for ff in time])
                elif time_format == 'BJD_UTC':
                    return np.array([self._bjd_tdb(None, curve_fit(self._bjd_utc, [0], [ff], p0=[ff])[0][0])
                                     for ff in time])
                elif time_format in ['HJD_TDB', 'HJD_TT']:
                    return np.array([self._bjd_tdb(None, curve_fit(self._hjd_tdb, [0], [ff], p0=[ff])[0][0])
                                     for ff in time])
                elif time_format == 'HJD_UTC':
                    return np.array([self._bjd_tdb(None, curve_fit(self._hjd_utc, [0], [ff], p0=[ff])[0][0])
                                     for ff in time])
                else:
                    raise PyLCInputError(
                        'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                        'HJD_UTC, HJD_BJD, HJD_TT')

            except:
                raise PyLCInputError('Not valid input for time')

    def convert_to_jd(self, time, time_format):
        """Convert times to JD_UTC.

        Parameters
        ----------
        time : float or array_like
            Time(s) to convert.
        time_format : str
            Time system of ``time``: ``'BJD_TDB'``, ``'BJD_TT'``,
            ``'BJD_UTC'``, ``'HJD_TDB'``, ``'HJD_TT'``, ``'HJD_UTC'``, ``'JD_UTC'`` or
            ``'MJD_UTC'``.

        Returns
        -------
        float or numpy.ndarray
            The time(s) in JD_UTC: a float for a float input, an array otherwise.

        Raises
        ------
        PyLCInputError
            If ``time_format`` is not valid or ``time`` cannot be converted to floats.

        """

        if isinstance(time, float):

            if time_format in ['BJD_TDB', 'BJD_TT']:
                return curve_fit(self._bjd_tdb, [0], [time], p0=[time])[0][0]
            elif time_format == 'JD_UTC':
                return time
            elif time_format == 'MJD_UTC':
                return time + 2400000.5
            elif time_format == 'BJD_UTC':
                return curve_fit(self._bjd_utc, [0], [time], p0=[time])[0][0]
            elif time_format in ['HJD_TDB', 'HJD_TT']:
                return curve_fit(self._hjd_tdb, [0], [time], p0=[time])[0][0]
            elif time_format == 'HJD_UTC':
                return curve_fit(self._hjd_utc, [0], [time], p0=[time])[0][0]
            else:
                raise PyLCInputError(
                    'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                    'HJD_UTC, HJD_BJD, HJD_TT')

        else:
            try:
                time = np.array(time, dtype=float)

                if time_format in ['BJD_TDB', 'BJD_TT']:
                    return np.array([curve_fit(self._bjd_tdb, [0], [ff], p0=[ff])[0][0] for ff in time])
                elif time_format == 'JD_UTC':
                    return time
                elif time_format == 'MJD_UTC':
                    return time + 2400000.5
                elif time_format == 'BJD_UTC':
                    return np.array([curve_fit(self._bjd_utc, [0], [ff], p0=[ff])[0][0] for ff in time])
                elif time_format in ['HJD_TDB', 'HJD_TT']:
                    return np.array([curve_fit(self._hjd_tdb, [0], [ff], p0=[ff])[0][0] for ff in time])
                elif time_format == 'HJD_UTC':
                    return np.array([curve_fit(self._hjd_utc, [0], [ff], p0=[ff])[0][0] for ff in time])
                else:
                    raise PyLCInputError(
                        'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                        'HJD_UTC, HJD_BJD, HJD_TT')

            except:
                raise PyLCInputError('Not valid input for time')

    def convert_to_mjd(self, time, time_format):
        """Convert times to MJD_UTC.

        Parameters
        ----------
        time : float or array_like
            Time(s) to convert.
        time_format : str
            Time system of ``time``: ``'BJD_TDB'``, ``'BJD_TT'``,
            ``'BJD_UTC'``, ``'HJD_TDB'``, ``'HJD_TT'``, ``'HJD_UTC'``, ``'JD_UTC'`` or
            ``'MJD_UTC'``.

        Returns
        -------
        float or numpy.ndarray
            The time(s) in MJD_UTC (JD_UTC - 2400000.5).

        Raises
        ------
        PyLCInputError
            If ``time_format`` is not valid or ``time`` cannot be converted to floats.

        """

        return self.convert_to_jd(time, time_format) - 2400000.5

    def _hjd_utc(self, x, jd):

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(jd)

        a = ssb_d / 60.0 / 24.0
        b = self.dec.sin() * np.sin(ssb_dec)
        c = self.dec.cos() * np.cos(ssb_dec) * np.cos(self.ra.rad() - ssb_ra)

        return jd - a * (b + c)

    def _hjd_tdb(self, x, jd):

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(jd)

        a = ssb_d / 60.0 / 24.0
        b = self.dec.sin() * np.sin(ssb_dec)
        c = self.dec.cos() * np.cos(ssb_dec) * np.cos(self.ra.rad() - ssb_ra)

        return jd - a * (b + c) + ssb_dt / 60.0 / 60.0 / 24.0

    def _bjd_utc(self, x, jd):

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.barycentre(jd)

        a = ssb_d / 60.0 / 24.0
        b = self.dec.sin() * np.sin(ssb_dec)
        c = self.dec.cos() * np.cos(ssb_dec) * np.cos(self.ra.rad() - ssb_ra)

        return jd - a * (b + c)

    def _bjd_tdb(self, x, jd):

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.barycentre(jd)

        a = ssb_d / 60.0 / 24.0
        b = self.dec.sin() * np.sin(ssb_dec)
        c = self.dec.cos() * np.cos(ssb_dec) * np.cos(self.ra.rad() - ssb_ra)

        return jd - a * (b + c) + ssb_dt / 60.0 / 60.0 / 24.0

    def distance_on_sphere(self, other):
        """Calculate the angular distance to another target.

        Parameters
        ----------
        other : FixedTarget
            The other target.

        Returns
        -------
        Rad
            Angular distance between the two targets.

        Raises
        ------
        PyLCInputError
            If ``other`` is not a :class:`FixedTarget`.

        """

        _request_target(other)

        d1, d2, dra = self.dec.rad(), other.dec.rad(), (self.ra - other.ra).rad()
        num = np.hypot(np.cos(d2) * np.sin(dra),
                       np.cos(d1) * np.sin(d2) - np.sin(d1) * np.cos(d2) * np.cos(dra))
        den = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(dra)
        return Rad(np.arctan2(num, den))

    def __str__(self):
        return 'plc.FixedTarget(RA(hms)/DEC(dms): {0})'.format(self.coord)

    def __repr__(self):
        return self.__str__()


def _is_target(item):
    if isinstance(item, FixedTarget):
        return True
    else:
        return False


def _request_target(item):
    if _is_target(item):
        pass
    else:
        raise PyLCInputError('A plc.Target object is required (plc.FixedTarget)')
