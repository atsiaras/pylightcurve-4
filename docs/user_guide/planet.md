# Planets

In PyLightcurve 4.0, all the calculations related to an exoplanet are available
through the {class}`~pylightcurve.Planet` class:

- limb-darkening coefficients and planet-to-star flux ratios ({doc}`filters`),
- transit and eclipse models, depths and durations, and the eclipse mid-time
  ({doc}`models`),
- the orbital position, the planet-star projected distance and the orbital phase
  ({doc}`models`),
- transit and eclipse fitting ({doc}`fitting`).

## Planets from the catalogue

The quickest way to get a planet is from the Exoplanet Characterisation Catalogue
(ECC), developed as part of the [ExoClock project](https://www.exoclock.space):

```python
planet = plc.get_planet('HD209458b')
```

The search ignores case, spaces and hyphens, so `'hd 209458 b'` also works. To list
all the planets in the catalogue:

```python
all_planets = plc.get_all_planets()
```

The catalogue entries of the planet and of its host star are kept in
`planet.all_data`.

To find the catalogue planet at a position on the sky, use
{func}`~pylightcurve.locate_planet`:

```python
planet = plc.locate_planet(plc.Hours('22:03:10.77'), plc.Degrees('+18:53:03.5'))
```

It searches within 1 arcminute by default, and raises
{class}`~pylightcurve.PyLCLibraryError` if no planet is found.

## Planets with your own parameters

You can also create a planet with your own parameters. For HD 209458 b:

```python
planet = plc.Planet(
    name='HD209458b',
    ra=330.795,                    # degrees, or plc.Hours('22:03:10.7729')
    dec=18.884,                    # degrees, or plc.Degrees('+18:53:03.548')
    stellar_logg=4.36,             # log(g), cgs
    stellar_temperature=6065.0,    # Kelvin
    stellar_metallicity=0.0,       # dex
    rp_over_rs=0.12086,            # Rp/R*
    period=3.5247486,              # days
    sma_over_rs=8.76,              # a/R*
    eccentricity=0.0,
    inclination=86.71,             # degrees, or a plc.Degrees object
    periastron=0.0,                # degrees, or a plc.Degrees object
    mid_time=2452826.62928,
    mid_time_format='BJD_TDB',     # see "Angles and times"
    ldc_method='claret',           # default; or 'linear', 'quad', 'sqrt'
    ldc_stellar_model='phoenix',   # default; or 'atlas'
    albedo=0.15,                   # default
    emissivity=1.0,                # default
)
```

Angles (`ra`, `dec`, `inclination`, `periastron`) can be given as numbers, in
degrees, or as angle objects such as `plc.Hours('22:03:10.7729')` (see
{doc}`angles_and_times`).

The mid-time can be given in any of the supported time systems
({doc}`angles_and_times`); it is converted to BJD_TDB and stored in
`planet.mid_time`. The mid-time of the eclipse that follows it is calculated from
the orbit and stored in `planet.eclipse_mid_time`, also in BJD_TDB:

```python
print(planet.mid_time, planet.eclipse_mid_time)
```

The other parameters are stored as attributes with the same names, for example
`planet.period` and `planet.sma_over_rs`.
