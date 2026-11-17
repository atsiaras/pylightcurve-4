# Angles and times

## Angles

PyLightcurve 4.0 has three angle classes, which differ only in how the angle is
defined:

- {class}`~pylightcurve.Degrees`: degrees, arcminutes and arcseconds,
- {class}`~pylightcurve.Hours`: hours, minutes and seconds (1 hour = 15 degrees),
  as for right ascensions,
- {class}`~pylightcurve.Rad`: radians.

Degrees and hours can be given as numbers or as strings:

```python
ra = plc.Hours('22:03:10.7729')
ra = plc.Hours(22, 3, 10.7729)

dec = plc.Degrees('+18:53:03.548')
dec = plc.Degrees(18, 53, 3.548)

angle = plc.Rad(0.5)
```

Whichever way an angle is defined, it can be expressed in any format:

```python
ra = plc.Hours('22:03:10.7729')

ra.deg()          # degrees, 0 to 360
ra.hours()        # hours, 0 to 24
ra.rad()          # radians
ra.dms()          # degrees, arcminutes and arcseconds, printed as 'dd:mm:ss.s'
ra.hms()          # hours, minutes and seconds, printed as 'hh:mm:ss.s'
ra.deg_coord()    # degrees, -90 to 90, as for a declination
ra.dms_coord()    # '+dd:mm:ss.s', -90 to 90, as for a declination
```

Angles are wrapped into the range 0 to 360 degrees, so `plc.Degrees(-10)` is the
same angle as `plc.Degrees(350)`; use {meth}`~pylightcurve.Degrees.deg_coord` to
get negative declinations. Angles can be added and subtracted, multiplied and
divided by numbers, and have {meth}`~pylightcurve.Degrees.sin`,
{meth}`~pylightcurve.Degrees.cos` and {meth}`~pylightcurve.Degrees.tan` methods.
{func}`~pylightcurve.arcsin`, {func}`~pylightcurve.arccos` and
{func}`~pylightcurve.arctan` return angle objects, and `plc.pi` is 180 degrees.

## Time systems

PyLightcurve 4.0 supports these time systems:

| Name | Time system |
|:--|:--|
| `JD_UTC` | Julian Date at the observer, in UTC |
| `MJD_UTC` | Modified Julian Date (JD - 2400000.5) at the observer, in UTC |
| `HJD_UTC`, `HJD_TDB` | Heliocentric Julian Date (time at the centre of the Sun), in UTC or TDB |
| `BJD_UTC`, `BJD_TDB` | Barycentric Julian Date (time at the barycentre of the Solar System), in UTC or TDB |

`HJD_TT` and `BJD_TT` are also accepted, and treated as `HJD_TDB` and `BJD_TDB`.
The observer is assumed to be at the centre of the Earth.

Internally, PyLightcurve uses BJD_TDB. The {class}`~pylightcurve.Planet` methods
accept times in any of these systems, and convert them to BJD_TDB using the
position of the host star.

## Converting times

The light-travel time to the Sun or to the barycentre depends on the direction of
the target, so time conversions need a {class}`~pylightcurve.FixedTarget`:

```python
ra = plc.Hours('22:03:10.7729')
dec = plc.Degrees('+18:53:03.548')
target = plc.FixedTarget(ra, dec)

mid_time_in_hjd_utc = 2458485.00380255
mid_time_in_bjd_tdb = target.convert_to_bjd_tdb(mid_time_in_hjd_utc, 'HJD_UTC')
mid_time_in_jd_utc = target.convert_to_jd(mid_time_in_hjd_utc, 'HJD_UTC')
mid_time_in_mjd_utc = target.convert_to_mjd(mid_time_in_hjd_utc, 'HJD_UTC')
```

The conversions accept a single float or an array. A planet has its own target,
`planet.target`:

```python
planet = plc.get_planet('HD209458b')
bjd_tdb = planet.target.convert_to_bjd_tdb(2458485.00380255, 'HJD_UTC')
```

The conversions use ephemeris tables that PyLightcurve downloads on first use
(see {doc}`/installation`), so they only work for dates covered by these tables.

The angular distance between two targets is given by
{meth}`~pylightcurve.FixedTarget.distance_on_sphere`:

```python
other = plc.FixedTarget(plc.Hours('22:03:00'), plc.Degrees('+19:00:00'))
print(target.distance_on_sphere(other).deg())
```
