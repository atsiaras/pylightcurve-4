
__all__ = ['get_planet', 'get_all_planets', 'locate_planet']

import numpy as np

from pylightcurve.errors import *
from pylightcurve.__databases__ import plc_data
from pylightcurve.models.exoplanet import Planet
from pylightcurve.spacetime.angles import Degrees, Hours, _request_angle
from pylightcurve.spacetime.targets import FixedTarget


def _flat_name(name):

    flat_name_list = [
        [' ', ''],
        ['-', ''],
        ['cancri', 'cnc'],
        ['hatp10', 'wasp11'],
        ['wasp40', 'hatp27'],
        ['wasp51', 'hatp30'],
        ['wasp86', 'kelt12'],
        ['kelt22', 'wasp173'],
    ]

    name = name.lower()

    for char in flat_name_list:
        name = name.replace(char[0], char[1])

    return name


def _search_by_planet(name):

    planets = plc_data.ecc()['planets']

    name_or = name

    if name in planets:
        return name

    else:
        for i in planets:
            if _flat_name(i) == _flat_name(name):
                return str(i)
            elif (_flat_name(i)[-1] == _flat_name(name)[-1] and _flat_name(i)[:-2] == _flat_name(name)[:-1] and
                  _flat_name(i)[-2] in ['a', 'b', 'n']):
                return str(i)

    raise PyLCInputError('No planet {0} found in the catalogue.'.format(name_or))


def get_all_planets():
    """List the planets in the Exoplanet Characterisation Catalogue.

    Returns
    -------
    list of str
        Names of all the planets that :func:`get_planet` can create.

    Examples
    --------
    >>> import pylightcurve as plc
    >>> names = plc.get_all_planets()

    """

    return list(plc_data.ecc()['planets'].keys())


def get_planet(name):
    """Create a :class:`Planet` from the Exoplanet Characterisation Catalogue.

    The stellar and planetary parameters are taken from the Exoplanet
    Characterisation Catalogue (ECC), developed as part of the
    `ExoClock project <https://www.exoclock.space>`_
    (`Kokori et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020arXiv201207478K/abstract>`_).

    Parameters
    ----------
    name : str
        Name of the planet. The search ignores case, spaces and hyphens, and accepts
        some alternative names, so ``'HD209458b'``, ``'hd 209458 b'`` and
        ``'HD-209458b'`` are equivalent. See :func:`get_all_planets` for the available
        names.

    Returns
    -------
    Planet
        The planet, with the default options of :class:`Planet`. The full catalogue
        entries are stored in its ``all_data`` attribute.

    Raises
    ------
    PyLCInputError
        If the planet is not in the catalogue.

    Examples
    --------
    >>> import pylightcurve as plc
    >>> planet = plc.get_planet('HD209458b')
    >>> planet.transit_duration('COUSINS_R')

    """

    name = _search_by_planet(name)

    planet_data = plc_data.ecc()['planets'][name]
    star_data = plc_data.ecc()['stars'][planet_data['star']]

    planet = Planet(
        name,
        Hours(star_data['ra']),
        Degrees(star_data['dec']),
        planet_data['logg'],
        planet_data['teff'],
        planet_data['meta'],
        planet_data['rp_over_rs'],
        planet_data['ephem_period'],
        planet_data['sma_over_rs'],
        planet_data['eccentricity'],
        planet_data['inclination'],
        planet_data['periastron'],
        planet_data['ephem_mid_time'],
        'BJD_TDB',
    )

    planet.all_data = {'planet': planet_data, 'star': star_data}

    return planet


def locate_planet(ra, dec, radius=Degrees(0, 1, 0)):
    """Find the catalogue planet closest to a position on the sky.

    Parameters
    ----------
    ra : float or Degrees or Hours or Rad
        Right ascension. A number is in degrees.
    dec : float or Degrees or Hours or Rad
        Declination. A number is in degrees.
    radius : float or Degrees or Hours or Rad, optional
        Search radius. A number is in degrees. Default is 1 arcminute.

    Returns
    -------
    Planet
        The planet closest to the position, created with :func:`get_planet`.

    Raises
    ------
    PyLCLibraryError
        If no planet is within ``radius``.

    Notes
    -----
    This function creates every planet in the catalogue, so it takes some time.

    Examples
    --------
    >>> import pylightcurve as plc
    >>> planet = plc.locate_planet(plc.Hours('22:03:10.77'), plc.Degrees('+18:53:03.5'))

    """

    if isinstance(ra, (int, float, np.integer, np.floating)):
        ra = Degrees(float(ra))
    else:
        _request_angle(ra)

    if isinstance(dec, (int, float, np.integer, np.floating)):
        dec = Degrees(float(dec))
    else:
        _request_angle(dec)

    if isinstance(radius, (int, float, np.integer, np.floating)):
        radius = Degrees(float(radius))
    else:
        _request_angle(radius)

    pointing = FixedTarget(ra, dec)

    test_planets = []

    for test_planet_name in get_all_planets():
        test_planet = get_planet(test_planet_name).target
        test_planets.append([pointing.distance_on_sphere(test_planet).deg(), test_planet_name])

    test_planets.sort()

    if test_planets[0][0] < radius.deg():
        return get_planet(test_planets[0][1])
    else:
        raise PyLCLibraryError('Planet could not be located')
