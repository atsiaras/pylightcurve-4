
__all__ = ['EmceeFitting', 'values_to_print']

import emcee
import numpy as np
import warnings

from pylightcurve.errors import *
from pylightcurve.processes.counter import Counter
from pylightcurve.processes.files import save_dict
from pylightcurve.plots.plots_fitting import plot_mcmc_corner, plot_mcmc_traces, plot_mcmc_fitting
from pylightcurve.analysis.distributions import one_d_distribution


#  emcee


class EmceeFitting:
    """Fit a model to data with MCMC, using emcee.

    The priors are uniform within ``[limits1, limits2]``. A parameter whose limits
    are ``NaN`` is kept fixed at its initial value. Call :meth:`run_mcmc` to fit; the
    results are then in the ``results`` attribute.

    Parameters
    ----------
    input_data_x : array_like
        Independent variable, passed to ``model`` unchanged.
    input_data_y : array_like
        Data to fit.
    input_data_y_unc : array_like
        Uncertainties of ``input_data_y``.
    model : callable
        Model function, called as ``model(input_data_x, *parameters)`` with *all*
        the parameters (fixed and free) in the order of ``initials``.
    initials : array_like
        Initial value of every parameter.
    limits1, limits2 : array_like
        Lower and upper prior limits of every parameter. Use ``numpy.nan`` for both to
        keep a parameter fixed.
    walkers : int
        Number of MCMC walkers.
    iterations : int
        Total number of MCMC samples, summed over all walkers. Each walker makes
        ``iterations // walkers`` steps.
    burn_in : int
        Number of samples, summed over all walkers, discarded from the start of the
        chain.
    data_x_name, data_y_name : str, optional
        Keys under which the data are stored in ``results['input_series']``. Default
        is ``'x'`` and ``'y'``.
    data_x_print_name, data_y_print_name : str, optional
        Axis labels for the plots. Default is ``'x'`` and ``'y'``.
    parameters_names : list of str, optional
        Parameter names, used as keys in ``results['parameters']``. Default is
        ``'p0'``, ``'p1'``, and so on.
    parameters_print_names : list of str, optional
        Parameter labels for the plots (LaTeX allowed, without the dollar signs).
        Default is ``'p0'``, ``'p1'``, and so on.
    counter : str, optional
        Label of the progress counter. Default is ``'MCMC'``.
    strech_prior : float, optional
        The walkers start close to the initial values, within
        ``1 / strech_prior`` of the distance to each limit. Default is 1000.0.
    function_to_call : callable, optional
        Function called without arguments before every evaluation of the posterior,
        for example to update a progress display.

    Attributes
    ----------
    results : dict
        The results, filled by :meth:`run_mcmc`. See the Notes.

    Notes
    -----
    After :meth:`run_mcmc`, ``results`` contains:

    - ``'mcmc'``: the ``iterations``, ``walkers`` and ``burn_in`` used.
    - ``'input_series'``: the input data.
    - ``'parameters'``: a dict per parameter, with its ``value``, the lower and upper
      uncertainties ``m_error`` and ``p_error``, the same values rounded as strings
      (``print_value``, ``print_m_error``, ``print_p_error``), the ``trace`` of the
      samples, its histogram (``trace_bins``, ``trace_counts``), the ``initial``
      value and the prior limits (``min_allowed``, ``max_allowed``). For fixed
      parameters only ``value`` is set.
    - ``'parameters_final'``: the list of final values of all the parameters.
    - ``'output_series'``: the best-fit ``model`` and the ``residuals``.
    - ``'statistics'``: residual statistics (``res_mean``, ``res_std``,
      ``res_rms``, ``res_chi_sqr``, ``res_red_chi_sqr``, ``res_autocorr``,
      ``res_max_autocorr``) and the correlation matrix of the fitted parameters
      (``corr_matrix``, ``corr_variables``).

    The value of each parameter is the peak of its posterior histogram, and the
    uncertainties enclose 68% of the samples on each side of the peak (see
    :func:`one_d_distribution`). Samples further than five median absolute
    deviations from the median, in any parameter, are not used.

    Examples
    --------
    >>> import numpy as np
    >>> import pylightcurve as plc
    >>> x = np.linspace(0, 10, 100)
    >>> y = 2.0 * x + 1.0 + np.random.normal(0, 0.5, 100)
    >>> def line(x, slope, intercept):
    ...     return slope * x + intercept
    >>> fit = plc.EmceeFitting(x, y, np.full(100, 0.5), line, [1.0, 0.0],
    ...                        [0.0, -5.0], [5.0, 5.0], 20, 20000, 5000,
    ...                        parameters_names=['slope', 'intercept'])
    >>> fit.run_mcmc()
    >>> fit.results['parameters']['slope']['value']

    """

    def __init__(self, input_data_x, input_data_y, input_data_y_unc,
                 model, initials, limits1, limits2, walkers, iterations, burn_in,
                 data_x_name='x', data_y_name='y', data_x_print_name='x', data_y_print_name='y',
                 parameters_names=None, parameters_print_names=None, counter='MCMC', strech_prior=1000.0,
                 function_to_call=None):

        self.input_data_x = input_data_x
        self.input_data_y = input_data_y
        self.input_data_y_unc = input_data_y_unc

        self.data_x_name = data_x_name
        self.data_y_name = data_y_name
        self.data_x_print_name = data_x_print_name
        self.data_y_print_name = data_y_print_name

        self.model = model

        self.initials = np.array(initials, dtype=float)
        self.limits1 = np.array(limits1, dtype=float)
        self.limits2 = np.array(limits2, dtype=float)

        self.walkers = int(walkers)

        self.iterations_per_walker = int(int(iterations) / walkers)
        self.iterations = self.iterations_per_walker * walkers

        self.burn_in = int(burn_in)

        self.names = ['p{0}'.format(ff) for ff in range(len(initials))]
        if parameters_names:
            self.names = parameters_names

        self.print_names = ['p{0}'.format(ff) for ff in range(len(initials))]
        if parameters_print_names:
            self.print_names = parameters_print_names

        self.counter = Counter(counter, self.iterations + self.walkers, 100)

        self.strech_prior = strech_prior

        self.results = {
            'mcmc': {
                'iterations': self.iterations,
                'walkers': self.walkers,
                'burn_in': self.burn_in,
            },
            'input_series': {},
            'parameters': {},
            'parameters_final': [],
            'output_series': {},
            'statistics': {}}

        self.fitted_parameters = []

        self.mcmc_run_complete = False

        self.function_to_call = function_to_call

    def run_mcmc(self):
        """Run the MCMC and fill the ``results`` attribute."""

        fitted_parameters_indices = np.where(~np.isnan(self.limits1 * self.limits2))[0]

        dimensions = len(fitted_parameters_indices)

        internal_limits1 = self.limits1[fitted_parameters_indices]
        internal_limits2 = self.limits2[fitted_parameters_indices]
        internal_initials = self.initials[fitted_parameters_indices]

        walkers_initial_positions = np.random.uniform(
            (internal_initials -
             (internal_initials - internal_limits1) / self.strech_prior)[:, None] * np.ones(self.walkers),
            (internal_initials +
             (internal_limits2 - internal_initials) / self.strech_prior)[:, None] * np.ones(self.walkers))
        walkers_initial_positions = np.swapaxes(walkers_initial_positions, 0, 1)

        def internal_model(theta):
            parameters = self.initials
            parameters[fitted_parameters_indices] = theta
            return self.model(self.input_data_x, *parameters)

        def likelihood(theta):
            if ((internal_limits1 < theta) & (theta < internal_limits2)).all():
                chi = (self.input_data_y - internal_model(theta)) / self.input_data_y_unc
                return -0.5 * (np.sum(chi * chi) +
                               np.sum(np.log(2.0 * np.pi * (self.input_data_y_unc * self.input_data_y_unc))))
            else:
                return -np.inf

        def prior(theta):
            if ((internal_limits1 < theta) & (theta < internal_limits2)).all():
                return 0.0
            return -np.inf

        # probability

        def probability_core(theta):
            return prior(theta) + likelihood(theta)

        probability_core_function = probability_core
        if self.function_to_call:
            def probability_core_function(theta):
                self.function_to_call()
                return probability_core(theta)

        def probability(theta):
            self.counter.update()
            return probability_core_function(theta)

        # run sampler

        sampler = emcee.EnsembleSampler(self.walkers, dimensions, probability)
        sampler.run_mcmc(walkers_initial_positions, int(self.iterations) // int(self.walkers))
        mcmc_results = sampler.get_chain(flat=True)

        self.results['input_series'][self.data_x_name] = self.input_data_x
        self.results['input_series'][self.data_y_name] = self.input_data_y
        self.results['input_series']['{0}_unc'.format(self.data_y_name)] = self.input_data_y_unc
        self.results['input_series_x'] = self.data_x_name
        self.results['input_series_x_print'] = self.data_x_print_name
        self.results['input_series_y'] = self.data_y_name
        self.results['input_series_y_print'] = self.data_y_print_name
        self.results['input_series_y_unc'] = '{0}_unc'.format(self.data_y_name)

        trace_to_analyse = 0
        vars_check = 0
        for var in range(len(self.names)):

            if not np.isnan(self.limits1[var]):

                trace = mcmc_results[:, np.where(fitted_parameters_indices == var)[0][0]][self.burn_in:]

                median = np.median(trace)
                mad = np.sqrt(np.median((trace - median) ** 2))

                trace_to_analyse += (trace > (median - 5 * mad)) * (trace < (median + 5 * mad))

                vars_check += 1

        trace_to_analyse = np.where(trace_to_analyse == vars_check)

        for var in range(len(self.names)):

            if np.isnan(self.limits1[var]):

                variable = {'name': self.names[var], 'print_name': self.print_names[var],
                            'initial': None, 'min_allowed': None, 'max_allowed': None,
                            'trace': None, 'trace_bins': None, 'trace_counts': None,
                            'value': self.initials[var], 'm_error': None, 'p_error': None,
                            'print_value': self.initials[var], 'print_m_error': '-', 'print_p_error': '-'}

            else:
                trace = mcmc_results[:, np.where(fitted_parameters_indices == var)[0][0]][self.burn_in:]

                trace = trace[trace_to_analyse]

                bins, counts, value, m_error, p_error = \
                    one_d_distribution(trace, confidence_interval=0.68)

                print_value, print_m_error, print_p_error = values_to_print(value, m_error, p_error)

                variable = {'name': self.names[var], 'print_name': self.print_names[var],
                            'initial': self.initials[var],
                            'min_allowed': self.limits1[var], 'max_allowed': self.limits2[var],
                            'trace': trace, 'trace_bins': bins, 'trace_counts': counts,
                            'value': value, 'm_error': m_error, 'p_error': p_error,
                            'print_value': print_value, 'print_m_error': print_m_error, 'print_p_error': print_p_error}

                self.fitted_parameters.append(self.names[var])

            self.results['parameters'][self.names[var]] = variable
            self.results['parameters_final'].append(variable['value'])

        self.results['output_series']['model'] = self.model(self.input_data_x, *self.results['parameters_final'])
        self.results['output_series']['residuals'] = self.input_data_y - self.results['output_series']['model']

        to_correlate = []
        for parameter in self.fitted_parameters:
            to_correlate.append(self.results['parameters'][parameter]['trace'])
        correlation_matrix = np.corrcoef(to_correlate)

        res_autocorr = np.correlate(self.results['output_series']['residuals'],
                                    self.results['output_series']['residuals'], mode='full')
        res_autocorr = res_autocorr[res_autocorr.size // 2:] / res_autocorr[res_autocorr.size // 2:][0]

        self.results['statistics']['res_autocorr'] = res_autocorr
        self.results['statistics']['res_max_autocorr'] = np.max(res_autocorr[1:])
        self.results['statistics']['res_mean'] = np.mean(self.results['output_series']['residuals'])
        self.results['statistics']['res_std'] = np.std(self.results['output_series']['residuals'])
        self.results['statistics']['res_rms'] = np.sqrt(np.mean(self.results['output_series']['residuals']**2))
        self.results['statistics']['res_chi_sqr'] = np.sum(
            (self.results['output_series']['residuals'] ** 2) / (self.input_data_y_unc ** 2))
        self.results['statistics']['res_red_chi_sqr'] = (
                self.results['statistics']['res_chi_sqr'] / (len(self.input_data_y_unc) - len(self.fitted_parameters)))
        self.results['statistics']['corr_matrix'] = correlation_matrix
        self.results['statistics']['corr_variables'] = ','.join(self.fitted_parameters)

        self.mcmc_run_complete = True

    def save_all(self, export_file):
        """Save the ``results`` dictionary to a pickle file.

        Parameters
        ----------
        export_file : str
            Path of the pickle file. Open it with :func:`open_dict`.

        Raises
        ------
        PyLCProcessError
            If :meth:`run_mcmc` has not been run.

        """

        if not self.mcmc_run_complete:
            raise PyLCProcessError('MCMC not completed')

        save_dict(self.results, export_file)

    def save_results(self, export_file):
        """Save the fitted values and the residual statistics to a text file.

        Parameters
        ----------
        export_file : str
            Path of the text file.

        Raises
        ------
        PyLCProcessError
            If :meth:`run_mcmc` has not been run.

        """

        if not self.mcmc_run_complete:
            raise PyLCProcessError('MCMC not completed')

        w = open(export_file, 'w')

        w.write('# variable\tresult\tuncertainty\n')

        for i in self.names:
            w.write('{0}\t{1}\t-{2} +{3}\n'.format(self.results['parameters'][i]['name'],
                                                   self.results['parameters'][i]['print_value'],
                                                   self.results['parameters'][i]['print_m_error'],
                                                   self.results['parameters'][i]['print_p_error']))

        w.write('\n#Residuals:\n')
        w.write('#Mean: {0}\n'.format(self.results['statistics']['res_mean']))
        w.write('#STD: {0}\n'.format(self.results['statistics']['res_std']))
        w.write('#RMS: {0}\n'.format(self.results['statistics']['res_rms']))
        w.write('#Max auto-correlation: {0}\n'.format(self.results['statistics']['res_max_autocorr']))
        w.write('#Chi squared: {0}\n'.format(self.results['statistics']['res_chi_sqr']))
        w.write('#Reduced chi squared: {0}\n'.format(self.results['statistics']['res_red_chi_sqr']))

        w.close()

    def plot_fitting(self, export_file):
        """Plot the data, the best-fit model and the residuals.

        Parameters
        ----------
        export_file : str
            Path of the image file. The format is set by the extension, for example
            ``.pdf`` or ``.png``.

        """

        plot_mcmc_fitting(self, export_file)

    def plot_corner(self, export_file):
        """Plot the posterior distributions and the correlations of the fitted parameters.

        Parameters
        ----------
        export_file : str
            Path of the image file.

        """

        plot_mcmc_corner(self, export_file)

    def plot_traces(self, export_file):
        """Plot the MCMC traces of the fitted parameters.

        Parameters
        ----------
        export_file : str
            Path of the image file.

        """

        plot_mcmc_traces(self, export_file)


# decimal points and rounding


def values_to_print(value, error_minus, error_plus):
    """Round a value and its uncertainties for printing.

    The number of decimal places keeps two significant digits of the smaller
    uncertainty, with at least two decimal places.

    Parameters
    ----------
    value : float
        The value.
    error_minus, error_plus : float
        The lower and upper uncertainties.

    Returns
    -------
    tuple of str
        The rounded value, lower uncertainty and upper uncertainty.

    Examples
    --------
    >>> import pylightcurve as plc
    >>> plc.values_to_print(0.123456, 0.0012, 0.0034)
    ('0.1235', '0.0012', '0.0034')

    """

    value = float(value)
    error_minus = float(error_minus)
    error_plus = float(error_plus)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

    if error_minus >= 1.0 or error_minus == 0.0:
        digit1 = 2
    else:
        str_error_minus = '{0:.{test}f}'.format(error_minus, test=10 + abs(int(np.log10(error_minus))))
        digit1 = np.where([ff not in ['0', '.'] for ff in str_error_minus])[0][0]

    if error_plus >= 1.0 or error_plus == 0.0:
        digit2 = 2
    else:
        str_error_plus = '{0:.{test}f}'.format(error_plus, test=10 + abs(int(np.log10(error_plus))))
        digit2 = np.where([ff not in ['0', '.'] for ff in str_error_plus])[0][0]

    width = max(2, digit1, digit2)

    print_value = '{0:.{width}f}'.format(round(value, width), width=width)
    print_m_error = '{0:.{width}f}'.format(round(error_minus, width), width=width)
    print_p_error = '{0:.{width}f}'.format(round(error_plus, width), width=width)

    return print_value, print_m_error, print_p_error
