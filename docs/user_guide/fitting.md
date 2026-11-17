# Fitting light curves

A {class}`~pylightcurve.Planet` can fit transit or eclipse light curves from
different sources (with different time systems, flux formats and exposure times),
different filters (multi-colour) and different epochs (multi-epoch), all together.
The fit uses MCMC, through [emcee](https://emcee.readthedocs.io).

## Adding observations

First add the observations to the planet with
{meth}`~pylightcurve.Planet.add_observation`. Add either only transit observations
or only eclipse observations.

In this example we simulate two transit observations of HD 209458 b, in two
filters and two time systems. In practice you would load your own light curves.

```python
planet = plc.get_planet('HD209458b')

# first observation
time = np.arange(planet.mid_time - 0.1, planet.mid_time + 0.1, 0.001)

transit = planet.transit_integrated(time, time_format='BJD_TDB', exp_time=120,
                                    time_stamp='mid', filter_name='COUSINS_R',
                                    max_sub_exp_time=1)
systematics = 1.2 * (1 + 0.013 * (time - time[0]) + 0.03 * ((time - time[0]) ** 2))
error = np.random.normal(0, 0.002, len(time))
flux = transit * systematics + error
flux_unc = np.ones_like(error) * np.std(error)

planet.add_observation(
    time=time,                 # times of the data points
    time_format='BJD_TDB',     # time system of `time`
    exp_time=120,              # exposure time, in seconds
    time_stamp='mid',          # `time` is the start, middle or end of each exposure
    flux=flux,                 # flux values (or magnitudes)
    flux_unc=flux_unc,         # flux uncertainties
    flux_format='flux',        # 'flux' or 'mag'
    filter_name='COUSINS_R',   # a built-in filter or one added with add_filter
)

# second observation
time = np.arange(planet.mid_time - 0.05, planet.mid_time + 0.15, 0.001)

transit = planet.transit_integrated(time, time_format='HJD_UTC', exp_time=30,
                                    time_stamp='mid', filter_name='TESS',
                                    max_sub_exp_time=1)
systematics = 3.6 * (1 - 0.02 * (time - time[0]) + 0.05 * ((time - time[0]) ** 2))
error = np.random.normal(0, 0.0005, len(time))
flux = transit * systematics + error
flux_unc = np.ones_like(error) * np.std(error)

planet.add_observation(time=time, time_format='HJD_UTC', exp_time=30, time_stamp='mid',
                       flux=flux, flux_unc=flux_unc, flux_format='flux',
                       filter_name='TESS')
```

Each observation is stored in `planet.observations`, with its times converted to
BJD_TDB at mid-exposure. Magnitudes are converted to flux relative to the first data
point.

## Transit fitting

To fit all the observations together:

```python
planet.transit_fitting('hd209458b_fit')
```

The argument is the folder where the results are saved. Each observation is
modelled as

$$
F(t) = N \left(1 + L\,\Delta t + Q\,\Delta t^2\right) \times \mathrm{transit}(t)
$$

where $\Delta t$ is the time since the first data point of the observation, in
days, and $N$, $L$ and $Q$ are fitted separately for each observation. The transit
model is integrated over the exposure time of each observation, and uses the
limb-darkening coefficients and Rp/R* of its filter.

The main options are:

| Option | Default | Meaning |
|:--|:--|:--|
| `detrending_order` | 2 | Order of the polynomial trend: 0, 1 or 2 |
| `iterations` | 130000 | Total number of MCMC samples, summed over all walkers |
| `walkers` | 200 | Number of MCMC walkers |
| `burn_in` | 30000 | Number of samples, summed over all walkers, discarded at the start |
| `fit_rp_over_rs` | `True` | Fit Rp/R* |
| `fit_individual_rp_over_rs` | `True` | Fit a separate Rp/R* in each filter |
| `fit_sma_over_rs` | `False` | Fit a/R* |
| `fit_inclination` | `False` | Fit the inclination |
| `fit_mid_time` | `True` | Fit the mid-time |
| `fit_individual_times` | `True` | Fit a separate mid-time for each epoch |
| `fit_period` | `False` | Fit the period (requires `fit_individual_times=False`) |
| `fit_ldc1` to `fit_ldc4` | `False` | Fit the limb-darkening coefficients, per filter |
| `max_sub_exp_time` | 10 | Sub-exposure length for the integrated models, in seconds |
| `precision` | 3 | Numerical precision of the models, 0 to 6 |

The prior limits are set with:

| Option | Default | Meaning |
|:--|:--|:--|
| `fit_rp_over_rs_limits` | `[0.5, 2.0]` | Factors of the initial Rp/R* |
| `fit_sma_over_rs_limits` | `[0.5, 2.0]` | Factors of the initial a/R* |
| `fit_inclination_limits` | `[70.0, 90.0]` | Inclination, in degrees |
| `fit_mid_time_limits` | `[-0.2, 0.2]` | Days, relative to the predicted mid-time |
| `fit_period_limits` | `[0.8, 1.2]` | Factors of the initial period |
| `fit_ldc_limits` | `[0.0, 1.0]` | Limb-darkening coefficients |

The eccentricity and the argument of periastron are always fixed.

:::{tip}
The default number of samples gives reliable results but can take a long time. For
a first test, use fewer, for example
`planet.transit_fitting('test_fit', iterations=20000, walkers=50, burn_in=5000)`.
:::

## Eclipse fitting

For eclipse observations, use {meth}`~pylightcurve.Planet.eclipse_fitting`:

```python
planet.eclipse_fitting('my_eclipse_fit')
```

The model and the options are the same as for transits, except that:

- the planet-to-star flux ratio is fitted instead of the limb-darkening
  coefficients: `fit_fp_over_fs` (default `True`), `fit_individual_fp_over_fs`
  (default `False`) and `fit_fp_over_fs_limits` (default `[0.001, 1000.0]`, factors
  of the initial value);
- by default only Fp/F* is fitted: `fit_rp_over_rs`, `fit_individual_rp_over_rs`,
  `fit_mid_time` and `fit_individual_times` are all `False`.

## Results

The results are saved in the output folder (the folder is created if needed, and
existing files with the same names are overwritten):

| File | Contents |
|:--|:--|
| `results.txt` | Fitted values and uncertainties, and residual statistics |
| `results.pickle` | All the results, as a dictionary (open with {func}`~pylightcurve.open_dict`) |
| `diagnostics_dataset_<n>.txt` | Residual statistics of observation `n`, before and after detrending |
| `lightcurves_dataset_<n>.pdf` | Data, model and residuals of observation `n` |
| `correlations.pdf` | Posterior distributions and correlations of the fitted parameters |
| `traces.pdf` | MCMC traces of the fitted parameters |

The fitted value of each parameter is the peak of its posterior distribution, and
the uncertainties enclose 68% of the samples on each side of the peak.

The parameters are named as follows (`<n>` is the observation number, starting
from 1):

| Name | Parameter |
|:--|:--|
| `N_<n>`, `L_<n>`, `Q_<n>` | Normalisation and trend coefficients of observation `n` |
| `LDC1_<filter>` to `LDC4_<filter>` | Limb-darkening coefficients (transits) |
| `rp_<filter>`, or `rp` | Rp/R*, per filter or shared |
| `fp_<filter>`, or `fp` | Fp/F*, per filter or shared (eclipses) |
| `P`, `a`, `e`, `i`, `w` | Period, a/R*, eccentricity, inclination, argument of periastron |
| `T_0` | Mid-time at the reference epoch |
| `T_mid_<epoch>` | Transit mid-time of each epoch, with `fit_individual_times` |
| `T_<epoch>` | Eclipse mid-time of each epoch, with `fit_individual_times` |

To read the results in Python:

```python
results = plc.open_dict('hd209458b_fit/results.pickle')

rp = results['parameters']['rp_COUSINS_R']
print(rp['value'], rp['m_error'], rp['p_error'])
```

`results` has the keys described in {class}`~pylightcurve.EmceeFitting`, and
also:

- `'settings'`: the fitting options;
- `'detrended_input_series'`, `'detrended_output_series'` and
  `'detrended_statistics'`: the data, model and residuals divided by the trend,
  and their statistics;
- `'observations'`: the observations, each with its own model, residuals and
  statistics.

After the fit, the same information per observation is also in
`planet.observations`.

## General-purpose fitting

The MCMC fitting used by the planet is available for any model, through
{class}`~pylightcurve.EmceeFitting`:

```python
x = np.linspace(0, 10, 100)
y = 2.0 * x + 1.0 + np.random.normal(0, 0.5, 100)
y_unc = np.full(100, 0.5)

def line(x, slope, intercept):
    return slope * x + intercept

fitting = plc.EmceeFitting(
    x, y, y_unc, line,
    [1.0, 0.0],      # initial values
    [0.0, -5.0],     # lower limits
    [5.0, 5.0],      # upper limits
    20,              # walkers
    20000,           # total samples
    5000,            # burn-in samples
    parameters_names=['slope', 'intercept'],
)
fitting.run_mcmc()

print(fitting.results['parameters']['slope']['value'])
fitting.save_results('line_fit.txt')
fitting.plot_corner('line_fit_correlations.pdf')
```

To keep a parameter fixed, set both of its limits to `np.nan`.
