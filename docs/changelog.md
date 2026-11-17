# What's new

## 4.0.5

- New documentation, with docstrings for the whole public API.
- Fixed the MCMC burn-in: with emcee 3, `burn_in` now removes the first samples of
  the chain, instead of samples spread over the whole chain.
- Fixed `eclipse_integrated`, which could show a spurious dip near the transit
  phase.
- Fixed the classification of observations as transits or eclipses.
- `eclipse_fitting` no longer changes `planet.eclipse_mid_time`.
- Fixed the per-observation reduced chi-squared in the fitting diagnostics.
- Fixed the `force_same_rp_over_rs` and `force_same_fp_over_fs` fitting settings,
  which stored the opposite value.
- Plots: the corner plot includes fitted parameters whose initial value is 0, and
  the trace plots use the upper uncertainty for the upper limits.
- Limb-darkening coefficients for the Spitzer/IRAC filters always use the ATLAS
  models, also for non-solar metallicity.
- Angles can be given as integers in `Planet` and `locate_planet`.
- Fixed `FixedTarget.convert_to_jd` and `convert_to_mjd` for arrays of times, and
  made `distance_on_sphere` accurate for small distances.
- Fixed the printing of fractional seconds in `dms()` and `hms()`.
- Fixed `find_single_star` for stars beyond the number of rows of the image.
- `one_d_distribution` and `two_d_distribution` give a clear error for identical
  values, and the confidence intervals include the central bin correctly.
- Downloads verify SSL certificates (using `certifi`) and no longer leave partial
  files.
- Database fixes: clearer errors when a database is not available, and fixes in
  updating and extracting the databases.
- Tested with Python 3.12.

## 4.0.4

- Fixed packaging and test issues.
- Fixed LaTeX strings in the plots.
- Fixed a bug caused by the removal of `matplotlib.cm.get_cmap` in recent
  versions of matplotlib.

## 4.0.3

- Fixed a bug caused by the removal of `np.float` in recent versions of numpy.

## 4.0.2

- Database updates are checked in memory.
- Fixed the database updating loop.

## 4.0.1

- Fixed packaging and test issues.

## 4.0.0

### The `Planet` class

All the calculations related to an exoplanet are available through the new
{class}`~pylightcurve.Planet` class, which keeps the stellar and orbital parameters
and the filter-dependent parameters together. See {doc}`user_guide/planet`.

### New catalogue

PyLightcurve 4.0 no longer uses the Open Exoplanet Catalogue (OEC), because of the
large number of mistakes in it and the absence of parameter updates. OEC has been
replaced by the Exoplanet Characterisation Catalogue, developed as part of the
[ExoClock project](https://www.exoclock.space) and dedicated to transiting
exoplanets. See {func}`~pylightcurve.get_planet`.

### Multi-epoch and multi-colour fitting

Light curves from different sources (with different flux and time formats and
exposure times), different filters (multi-colour) and different epochs
(multi-epoch) can be fitted together. See {doc}`user_guide/fitting`.

### Eclipses

Eclipse models and eclipse fitting, and tools to calculate the eclipse mid-time,
depth and duration from the orbit of the planet. See {doc}`user_guide/models`.

### Limb-darkening coefficients

The default limb-darkening coefficients are calculated with
[ExoTETHyS](https://github.com/ucl-exoplanets/ExoTETHyS), which also covers stars
cooler than 3500 K. See {doc}`user_guide/filters`.

### Angles and times

New objects to convert between angle formats (for example `'+47:12:34.05'` to
degrees) and between time systems (for example HJD_UTC to BJD_TDB). See
{doc}`user_guide/angles_and_times`.
