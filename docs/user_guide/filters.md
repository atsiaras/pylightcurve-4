# Filters and limb darkening

Some parameters of a planet depend on the observing filter. They are needed to
calculate transit and eclipse models at different wavelengths, and to fit
light curves from different filters together:

- the planet-to-star radius ratio, `rp_over_rs`,
- the limb-darkening coefficients, `limb_darkening_coefficients`,
- the planet-to-star flux ratio, `fp_over_fs`.

## Built-in filters

PyLightcurve 4.0 includes these filters:

| Filters | Names |
|:--|:--|
| Unfiltered | `clear`, `luminance` |
| Johnson-Cousins | `JOHNSON_U`, `JOHNSON_B`, `JOHNSON_V`, `COUSINS_R`, `COUSINS_I` |
| 2MASS | `2mass_j`, `2mass_h`, `2mass_ks` |
| SDSS | `sdss_u`, `sdss_g`, `sdss_r`, `sdss_i`, `sdss_z` |
| Space telescopes | `Kepler`, `TESS` |
| Spitzer/IRAC | `irac1`, `irac2`, `irac3`, `irac4` |

The filter names are case-sensitive. To list the filters in your installation:

```python
from pylightcurve.__databases__ import plc_data
print(plc_data.all_filters())
```

## Default values

The first time a filter is used with a planet, PyLightcurve calculates its
parameters:

- `rp_over_rs` is the value given when the planet was created.
- The limb-darkening coefficients are interpolated from grids calculated with
  [ExoTETHyS](https://github.com/ucl-exoplanets/ExoTETHyS), for the stellar
  parameters, the `ldc_method` and the `ldc_stellar_model` of the planet. See
  {func}`~pylightcurve.exotethys`.
- `fp_over_fs` is the sum of the reflected light and the thermal emission of the
  planet, from its `albedo` and `emissivity`, assuming that the star and the
  planet emit as black bodies. See {func}`~pylightcurve.fp_over_fs`.

You can see these values with {meth}`Planet.filter <pylightcurve.Planet.filter>`:

```python
planet = plc.get_planet('HD209458b')

limb_darkening_coefficients = planet.filter('COUSINS_R').limb_darkening_coefficients
fp_over_fs = planet.filter('COUSINS_R').fp_over_fs
rp_over_rs = planet.filter('COUSINS_R').rp_over_rs   # same as planet.rp_over_rs
```

## Your own values

To define a new filter, or to replace the default values of a built-in one, use
{meth}`Planet.add_filter <pylightcurve.Planet.add_filter>`:

```python
planet.add_filter('my_filter', rp_over_rs, ldc1, ldc2, ldc3, ldc4, fp_over_fs)
```

The new name can then be used in every method that takes a filter name. The
limb-darkening coefficients must be for the `ldc_method` of the planet; give four
values even for laws with fewer coefficients, padding with zeros.

## Limb-darkening laws

The limb-darkening law is chosen with the `ldc_method` argument of
{class}`~pylightcurve.Planet`. With $\mu$ the cosine of the angle between the line
of sight and the normal to the stellar surface:

| `ldc_method` | Law |
|:--|:--|
| `'claret'` (default) | $I(\mu) = 1 - \sum_{k=1}^{4} a_k \left(1 - \mu^{k/2}\right)$ |
| `'quad'` | $I(\mu) = 1 - a_1 (1 - \mu) - a_2 (1 - \mu)^2$ |
| `'sqrt'` | $I(\mu) = 1 - a_1 (1 - \mu) - a_2 \left(1 - \sqrt{\mu}\right)$ |
| `'linear'` | $I(\mu) = 1 - a_1 (1 - \mu)$ |

## Stellar models

The limb-darkening grids come from two sets of stellar models, chosen with the
`ldc_stellar_model` argument of {class}`~pylightcurve.Planet`:

- `'phoenix'` (default): PHOENIX models, which also cover stars cooler than 3500 K
  but only solar metallicity. The metallicity is set to 0, with a message.
- `'atlas'`: ATLAS models, from 3500 K.

For the Spitzer/IRAC filters only the ATLAS grid is available, so PyLightcurve
switches to `'atlas'` automatically, with a message.
