# PyLightcurve 4.0

**Latest version v4.0.5**

**A Python package for modelling and analysing exoplanet light curves.**

```{image} _static/logo.png
:alt: PyLightcurve logo
:width: 200px
:align: right
```

:::{note}
This is the documentation of **PyLightcurve 4.0**, the long-term support version.
New projects should use the latest version, on the 
<a href="https://github.com/ucl-exoplanets/pylightcurve" target="_blank">main Pylightcurve repository</a>
:::

PyLightcurve 4.0 provides tools for:

- **Transit and eclipse models**: limb-darkened transit models, eclipse models and
  models integrated over the exposure time.
- **Light-curve fitting**: MCMC fitting of multi-epoch and multi-colour transit and
  eclipse light curves, with a polynomial trend for each observation.
- **Planet properties**: transit and eclipse depths, durations and mid-times, and
  orbital positions and phases.
- **Limb darkening**: limb-darkening coefficients for the built-in filters, from
  grids calculated with [ExoTETHyS](https://github.com/ucl-exoplanets/ExoTETHyS).
- **Angles and times**: conversions between angle formats, and between the JD,
  MJD, HJD and BJD time systems in UTC and TDB.

PyLightcurve is developed by [Angelos Tsiaras](https://www.angelostsiaras.com) and
collaborators.

```python
# Plot your first transit model!
import numpy as np
import pylightcurve as plc
import matplotlib.pyplot as plt

planet = plc.get_planet('HD209458b')
time = np.arange(planet.mid_time - 0.2, planet.mid_time + 0.2, 0.001)
model = planet.transit(time, 'BJD_TDB', 'COUSINS_R')

plt.figure(figsize=(9,4))
plt.plot((time-planet.mid_time)*24, model, 'k-', lw=2)
plt.ylabel('Relative flux')
plt.xlabel('Time [hours]')
```
```{image} _static/transit_model.png
:alt: transit model
:align: center
```

::::{grid} 2
:gutter: 3

:::{grid-item-card} Installation
:link: installation
:link-type: doc

Install PyLightcurve.
:::

:::{grid-item-card} Getting started
:link: quickstart
:link-type: doc

Compute your first transit model.
:::

:::{grid-item-card} User guide
:link: user_guide/index
:link-type: doc

The `Planet` class, filters, models, fitting, and angles and times.
:::

:::{grid-item-card} API reference
:link: api/index
:link-type: doc

Every public class and function.
:::

:::{grid-item-card} What's new in 4.0
:link: changelog
:link-type: doc

New features in 4.0 and the changes in each 4.0.x release.
:::
::::

```{toctree}
:maxdepth: 2
:hidden:

installation
quickstart
user_guide/index
api/index
changelog
citing
```
