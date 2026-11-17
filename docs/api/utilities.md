# Utilities

Helper functions that PyLightcurve uses internally. They are part of the public
package, and you can use them in your own analysis.

```{eval-rst}
.. currentmodule:: pylightcurve
```

## Distributions

```{eval-rst}
.. autofunction:: one_d_distribution
.. autofunction:: two_d_distribution
```

## Gaussian models and fitting

```{eval-rst}
.. autofunction:: gaussian
.. autofunction:: fit_gaussian
.. autofunction:: two_d_gaussian
.. autofunction:: fit_two_d_gaussian
.. autofunction:: curve_fit
```

## Numerical integration

```{eval-rst}
.. autofunction:: gauss_numerical_integration
.. autofunction:: sample_function
```

## Images

```{eval-rst}
.. autofunction:: find_single_star
```

## Files

```{eval-rst}
.. autofunction:: open_dict
.. autofunction:: save_dict
.. autofunction:: copy_dict
.. autofunction:: open_yaml
.. autofunction:: save_yaml
.. autofunction:: copy_yaml
.. autofunction:: open_fits
.. autofunction:: save_fits
.. autofunction:: copy_fits
.. autofunction:: fits_sci
.. autofunction:: fits_err
.. autofunction:: fits_sci_err
.. autofunction:: download
```

## Progress counter

```{eval-rst}
.. autoclass:: Counter
   :members: update
```
