# Quickstart

This page shows the main features of PyLightcurve 4.0 in a few lines. The
{doc}`user_guide/index` explains each of them in more detail.

```python
import numpy as np
import matplotlib.pyplot as plt
import pylightcurve as plc
```

## Get a planet

Create a {class}`~pylightcurve.Planet` from the Exoplanet Characterisation
Catalogue:

```python
planet = plc.get_planet('HD209458b')

print(planet.period, planet.rp_over_rs, planet.sma_over_rs, planet.inclination)
print(planet.mid_time)             # transit mid-time, BJD_TDB
print(planet.eclipse_mid_time)     # eclipse mid-time, BJD_TDB
```

## Model a transit and an eclipse

Every model needs the times, their time system and a filter:

```python
time = np.arange(planet.mid_time - 0.15, planet.mid_time + 0.15, 0.001)
transit = planet.transit(time, 'BJD_TDB', 'COUSINS_R')

eclipse_time = np.arange(planet.eclipse_mid_time - 0.15, planet.eclipse_mid_time + 0.15, 0.001)
eclipse = planet.eclipse(eclipse_time, 'BJD_TDB', 'irac1')

plt.plot(time - planet.mid_time, transit)
plt.xlabel('time from mid-transit [days]')
plt.ylabel('relative flux')
plt.show()
```

For observed data, use the models integrated over the exposure time:

```python
transit = planet.transit_integrated(time, 'BJD_TDB', exp_time=120, time_stamp='mid',
                                    filter_name='COUSINS_R')
```

## Transit and eclipse properties

```python
print(planet.transit_duration('COUSINS_R'))   # days
print(planet.transit_depth('COUSINS_R'))
print(planet.eclipse_duration('irac1'))       # days
print(planet.eclipse_depth('irac1'))
```

## Fit a light curve

Add one or more observations and fit them together:

```python
planet.add_observation(time=time, time_format='BJD_TDB', exp_time=120, time_stamp='mid',
                       flux=flux, flux_unc=flux_unc, flux_format='flux',
                       filter_name='COUSINS_R')

planet.transit_fitting('my_fit')
```

The results are saved in the folder `my_fit`. See {doc}`user_guide/fitting`.

## Convert times

```python
target = plc.FixedTarget(plc.Hours('22:03:10.7729'), plc.Degrees('+18:53:03.548'))
bjd_tdb = target.convert_to_bjd_tdb(2458485.00380255, 'HJD_UTC')
```
