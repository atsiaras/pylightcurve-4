# Model functions

Functions that take all the stellar and orbital parameters directly. The
{class}`~pylightcurve.Planet` methods call them with the parameters of the planet.
Times are in BJD_TDB.

```{eval-rst}
.. currentmodule:: pylightcurve
```

## Orbit

```{eval-rst}
.. autofunction:: planet_orbit
.. autofunction:: planet_star_projected_distance
.. autofunction:: planet_phase
```

## Transits

```{eval-rst}
.. autofunction:: transit
.. autofunction:: transit_integrated
.. autofunction:: transit_depth
.. autofunction:: transit_duration
```

## Eclipses

```{eval-rst}
.. autofunction:: eclipse
.. autofunction:: eclipse_integrated
.. autofunction:: eclipse_depth
.. autofunction:: eclipse_duration
.. autofunction:: eclipse_mid_time
```

## Limb darkening and planet emission

```{eval-rst}
.. autofunction:: exotethys
.. autofunction:: fp_over_fs
```
