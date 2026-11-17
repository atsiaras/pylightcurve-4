# Transit and eclipse models

## Models from a planet

For a series of times in any of the supported time systems
({doc}`angles_and_times`), a {class}`~pylightcurve.Planet` calculates transit and
eclipse models in any filter ({doc}`filters`):

```python
planet = plc.get_planet('HD209458b')

time_array = np.arange(planet.mid_time - 0.1, planet.mid_time + 0.1, 0.001)
transit_model = planet.transit(time_array, time_format='BJD_TDB', filter_name='COUSINS_R')

time_array = np.arange(planet.eclipse_mid_time - 0.1, planet.eclipse_mid_time + 0.1, 0.001)
eclipse_model = planet.eclipse(time_array, time_format='BJD_TDB', filter_name='COUSINS_R')
```

The transit model is the relative flux of the star, and the eclipse model the
relative flux of the system: both are 1 outside the transit or eclipse. The eclipse
model treats the planet as a uniform disk.

## Models integrated over the exposure time

The models above give the flux at each time. Observed data are integrated over the
exposure time, which smooths the light curve, especially for long exposures. For
these, use the integrated models, with the exposure time in seconds and the
position of the time stamps within each exposure (`'start'`, `'mid'` or `'end'`):

```python
transit_model_integrated = planet.transit_integrated(
    time_array, time_format='BJD_TDB', exp_time=120, time_stamp='mid',
    filter_name='COUSINS_R')

eclipse_model_integrated = planet.eclipse_integrated(
    time_array, time_format='BJD_TDB', exp_time=120, time_stamp='mid',
    filter_name='COUSINS_R')
```

Each exposure is split into sub-exposures, the model is calculated at the middle of
each, and the results are averaged. The length of the sub-exposures is set by
`max_sub_exp_time` (in seconds, default 10). Shorter sub-exposures give a more
precise model but take longer to calculate; the default is sufficient for all
known exoplanets.

## Numerical precision

The transit model is calculated by numerical integration over the stellar disk. The
`precision` argument, from 0 to 6, sets the number of integration points (2 for 0,
and 10 times `precision` otherwise). Higher values are more precise but slower. The
default for the models is 3, and for the depths 6.

## Transit and eclipse properties

The transit and eclipse depths and durations depend on the filter, through
`rp_over_rs`, the limb-darkening coefficients and `fp_over_fs`:

```python
transit_duration = planet.transit_duration('COUSINS_R')   # days
transit_depth = planet.transit_depth('COUSINS_R')
eclipse_duration = planet.eclipse_duration('COUSINS_R')   # days
eclipse_depth = planet.eclipse_depth('COUSINS_R')
```

The durations are from the first to the fourth contact (T14). The transit depth
includes limb darkening, so it is generally not equal to (Rp/R*)².

The eclipse mid-time, in BJD_TDB, is calculated when the planet is created:

```python
eclipse_mid_time = planet.eclipse_mid_time
```

## Orbit

For a series of times, a planet also calculates:

1. Its position on the orbit, with {meth}`~pylightcurve.Planet.planet_orbit`, as
   three arrays with the x, y and z coordinates. The star is at (0, 0, 0) and the
   observer at (+∞, 0, 0), and the planet is at periastron at mid-transit when the
   argument of periastron is 90 degrees. The coordinates are in units of the
   stellar radius.
2. The planet-star distance projected on the sky, with
   {meth}`~pylightcurve.Planet.planet_star_projected_distance`, in units of the
   stellar radius.
3. The orbital phase, with {meth}`~pylightcurve.Planet.planet_phase`: 0 at
   mid-transit, increasing by 1 every orbit.

```python
time_array = np.arange(planet.mid_time - 0.1, planet.mid_time + 0.1, 0.001)

x, y, z = planet.planet_orbit(time_array, 'BJD_TDB')
projected_distance = planet.planet_star_projected_distance(time_array, 'BJD_TDB')
phase = planet.planet_phase(time_array, 'BJD_TDB')
```

## Models without a planet

The same calculations are available as functions that take all the parameters
directly. They are faster when you need many models with different parameters,
for example in your own fitting code. Times are in BJD_TDB and the
limb-darkening coefficients are always given as four values:

```python
time_array = np.arange(2452826.52928, 2452826.72928, 0.001)

transit_model = plc.transit(
    [0.6, -0.4, 1.0, -0.4],    # limb-darkening coefficients
    0.12086,                   # rp_over_rs
    3.5247486,                 # period, days
    8.76,                      # sma_over_rs
    0.0,                       # eccentricity
    86.71,                     # inclination, degrees
    0.0,                       # argument of periastron, degrees
    2452826.62928,             # mid-time, BJD_TDB
    time_array,
    method='claret',
)
```

The functions are {func}`~pylightcurve.transit`,
{func}`~pylightcurve.transit_integrated`, {func}`~pylightcurve.transit_depth`,
{func}`~pylightcurve.transit_duration`, {func}`~pylightcurve.eclipse`,
{func}`~pylightcurve.eclipse_integrated`, {func}`~pylightcurve.eclipse_depth`,
{func}`~pylightcurve.eclipse_duration`, {func}`~pylightcurve.eclipse_mid_time`,
{func}`~pylightcurve.planet_orbit`,
{func}`~pylightcurve.planet_star_projected_distance`,
{func}`~pylightcurve.planet_phase`, {func}`~pylightcurve.exotethys` and
{func}`~pylightcurve.fp_over_fs`. See {doc}`/api/models`.

:::{note}
The integrated functions, {func}`~pylightcurve.transit_integrated` and
{func}`~pylightcurve.eclipse_integrated`, expect mid-exposure times. The
`time_stamp` option is only available in the {class}`~pylightcurve.Planet` methods.
:::
