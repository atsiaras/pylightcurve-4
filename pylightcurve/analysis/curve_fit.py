
__all__ = ['curve_fit']

import warnings

from scipy.optimize import curve_fit as scipy_curve_fit


def curve_fit(*args, **kwargs):
    """Wrap :func:`scipy.optimize.curve_fit`, hiding the covariance warning.

    All arguments are passed to :func:`scipy.optimize.curve_fit` unchanged. The only
    difference is that the "Covariance of the parameters could not be estimated"
    warning is not shown.

    Returns
    -------
    popt, pcov
        As returned by :func:`scipy.optimize.curve_fit`.

    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='Covariance of the parameters could not be estimated')
        return scipy_curve_fit(*args, **kwargs)
