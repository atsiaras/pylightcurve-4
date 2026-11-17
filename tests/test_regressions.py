"""Regression tests for the bugs fixed in 4.0.5 (logic review of 2026-11-16).

Each test fails on the 4.0.4 code and passes on the fixed code. The number in
each test name refers to the item in the review.

Most tests run offline. The ones that convert times, read the catalogue or
interpolate limb-darkening coefficients need the pylightcurve databases, like
the rest of the test suite.
"""

import os

import numpy as np
import pytest

import pylightcurve as plc
from pylightcurve.models.exoplanet_lc import eclipse, eclipse_integrated, eclipse_mid_time, transit_integrated, \
    exotethys


# HD 209458 b, rounded
RA, DEC = 330.795, 18.884
LOGG, TEFF, METALLICITY = 4.36, 6065.0, 0.0
RP, PERIOD, SMA, ECC, INC, W = 0.12086, 3.52474859, 8.76, 0.0, 86.71, 0.0
MID_TIME = 2452826.628521
LDCS = [0.4, 0.2, 0.1, 0.05]
FP = 0.0002


def make_planet():
    planet = plc.Planet('HD209458b', RA, DEC, LOGG, TEFF, METALLICITY,
                        RP, PERIOD, SMA, ECC, INC, W, MID_TIME, 'BJD_TDB')
    planet.add_filter('test', RP, *LDCS, FP)
    return planet


def line(x, slope, constant):
    return slope * x + constant


# ---------------------------------------------------------------------------
# 0, 0b: MCMC burn-in with emcee 3, integer initial values
# ---------------------------------------------------------------------------

class _StepLabelledSampler:
    """Stands in for emcee.EnsembleSampler. Every sample is labelled with its step number, and
    get_chain(flat=True) returns the samples in emcee 3's documented order: all the walkers of
    step 0, then all the walkers of step 1, and so on."""

    def __init__(self, nwalkers, ndim, log_prob_fn):
        self.nwalkers, self.ndim = nwalkers, ndim

    def run_mcmc(self, initial_state, nsteps):
        steps = np.arange(nsteps, dtype=float)[:, None, None]
        walkers = np.arange(self.nwalkers, dtype=float)[None, :, None] / 1000.0
        self.chain = np.broadcast_to(steps + walkers, (nsteps, self.nwalkers, self.ndim)).copy()

    def get_chain(self, flat=False):
        return self.chain.reshape(-1, self.ndim) if flat else self.chain


def test_00_burn_in_removes_the_first_steps(monkeypatch):
    import pylightcurve.analysis.optimisation as optimisation
    monkeypatch.setattr(optimisation.emcee, 'EnsembleSampler', _StepLabelledSampler)

    walkers, steps, burn_in_steps = 20, 100, 30
    x = np.arange(10.0)
    fit = plc.EmceeFitting(x, x, np.ones_like(x), line, [0.0, 0.0], [-1000, -1000], [1000, 1000],
                           walkers, walkers * steps, walkers * burn_in_steps, parameters_names=['s', 'c'])
    fit.run_mcmc()

    for name in ['s', 'c']:
        trace = fit.results['parameters'][name]['trace']
        assert len(trace) == walkers * (steps - burn_in_steps)
        assert np.floor(trace).min() == burn_in_steps
        assert np.floor(trace).max() == steps - 1


@pytest.mark.parametrize('initials', [[4.0, 0.0], [4, 0]])
def test_00_emcee_fit_recovers_least_squares_solution(initials):
    rng = np.random.default_rng(0)
    x = np.arange(0, 10, 0.1)
    unc = np.full_like(x, 0.5)
    y = line(x, 1.3, 2.3) + rng.normal(0, unc)

    design = np.vstack([x, np.ones_like(x)]).T
    exact = np.linalg.solve(design.T @ design, design.T @ y)
    exact_errors = np.sqrt(np.diag(np.linalg.inv(design.T @ design / unc[0] ** 2)))

    walkers = 50
    fit = plc.EmceeFitting(x, y, unc, line, initials, [-10, -10], [10, 10],
                           walkers, walkers * 600, walkers * 200, parameters_names=['s', 'c'])
    fit.run_mcmc()

    for name, value, error in zip(['s', 'c'], exact, exact_errors):
        result = fit.results['parameters'][name]
        assert abs(result['value'] - value) < 2 * error
        assert 0.7 < result['m_error'] / error < 1.3
        assert 0.7 < result['p_error'] / error < 1.3


# ---------------------------------------------------------------------------
# 1, 12: one_d_distribution
# ---------------------------------------------------------------------------

def test_01_confidence_interval_is_symmetric_for_symmetric_distribution():
    m_errors, p_errors = [], []
    for seed in range(20):
        sample = np.random.default_rng(seed).normal(0, 1, 100000)
        m_error, p_error = plc.one_d_distribution(sample, confidence_interval=0.68)[3:]
        m_errors.append(m_error)
        p_errors.append(p_error)

    # 4.0.4 gave 0.91 and 1.09
    assert 0.95 < np.mean(m_errors) < 1.05
    assert 0.95 < np.mean(p_errors) < 1.05
    assert 0.95 < np.mean(m_errors) / np.mean(p_errors) < 1.05


@pytest.mark.filterwarnings('ignore::RuntimeWarning')
def test_12_distribution_of_identical_values_gives_clear_error():
    with pytest.raises(plc.PyLCProcessError):
        plc.one_d_distribution(np.ones(50), confidence_interval=0.68)


def test_12_distribution_with_mostly_identical_values_does_not_crash():
    plc.one_d_distribution(np.array([1.0, 1, 1, 1, 2, 3]), confidence_interval=0.68)


# ---------------------------------------------------------------------------
# 2: eclipse_integrated
# ---------------------------------------------------------------------------

def test_02_eclipse_integrated_has_no_dip_outside_eclipse():
    emt = eclipse_mid_time(PERIOD, SMA, ECC, INC, W, MID_TIME)
    for time in [MID_TIME, MID_TIME + PERIOD / 4, MID_TIME - PERIOD / 4]:
        flux = eclipse_integrated(FP, RP, PERIOD, SMA, ECC, INC, W, emt, np.array([time]), 30)[0]
        assert flux == pytest.approx(1.0, abs=1e-12)


def test_02_eclipse_integrated_matches_eclipse_at_mid_eclipse():
    emt = eclipse_mid_time(PERIOD, SMA, ECC, INC, W, MID_TIME)
    time = np.array([emt])
    assert eclipse_integrated(FP, RP, PERIOD, SMA, ECC, INC, W, emt, time, 1)[0] == \
        pytest.approx(eclipse(FP, RP, PERIOD, SMA, ECC, INC, W, emt, time)[0], abs=1e-9)


# ---------------------------------------------------------------------------
# 4, 9, 15: eclipse_fitting (needs the databases)
# ---------------------------------------------------------------------------

def test_04_09_15_eclipse_fitting(tmp_path):
    planet = make_planet()
    original_eclipse_mid_time = planet.eclipse_mid_time

    rng = np.random.default_rng(1)
    centre = planet.eclipse_mid_time + 5 * PERIOD
    time = np.arange(centre - 0.12, centre + 0.12, 0.002)
    flux = planet.eclipse_integrated(time, 'BJD_TDB', 60, 'mid', 'test')
    flux_unc = np.full_like(flux, 0.0001)
    planet.add_observation(time, 'BJD_TDB', 60, 'mid', flux + rng.normal(0, flux_unc), flux_unc, 'flux', 'test')

    output_folder = str(tmp_path / 'eclipse')
    planet.eclipse_fitting(output_folder, iterations=20000, walkers=20, burn_in=5000)
    results = plc.open_dict(os.path.join(output_folder, 'results.pickle'))

    # 4: the planet is not changed by the fit
    assert planet.eclipse_mid_time == original_eclipse_mid_time

    # 9: reduced chi squared uses the number of fitted parameters
    fitted = [name for name in results['parameters'] if results['parameters'][name]['trace'] is not None]
    observation = planet.observations[0]
    assert observation['res_red_chi_sqr'] == pytest.approx(
        observation['res_chi_sqr'] / (len(observation['flux']) - len(fitted)))

    # 15: settings stored the right way round (defaults: one fp and one rp for all filters)
    assert results['settings']['force_same_fp_over_fs'] is True
    assert results['settings']['force_same_rp_over_rs'] is True


# ---------------------------------------------------------------------------
# 5: convert_to_jd with arrays (needs the databases)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('time_format', ['BJD_TDB', 'BJD_UTC', 'HJD_TDB', 'HJD_UTC'])
def test_05_convert_to_jd_accepts_arrays(time_format):
    target = plc.FixedTarget(plc.Degrees(RA), plc.Degrees(DEC))
    times = [2458485.0046340, 2458485.0146340, 2458486.5]

    converted = target.convert_to_jd(times, time_format)

    assert converted == pytest.approx([target.convert_to_jd(t, time_format) for t in times], abs=1e-9)


# ---------------------------------------------------------------------------
# 6: add_observation (needs the databases)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('offset_minutes', [-10, 10])
@pytest.mark.parametrize('kind', ['transit', 'eclipse'])
def test_06_observation_type_and_epoch(kind, offset_minutes):
    planet = make_planet()
    reference = planet.mid_time if kind == 'transit' else planet.eclipse_mid_time
    centre = reference + 100 * PERIOD + offset_minutes / 1440.0
    time = np.arange(centre - 0.12, centre + 0.12, 0.001)
    if kind == 'transit':
        flux = transit_integrated(LDCS, RP, PERIOD, SMA, ECC, INC, W, planet.mid_time, time, 30)
    else:
        flux = eclipse_integrated(FP, RP, PERIOD, SMA, ECC, INC, W, planet.eclipse_mid_time, time, 30)

    planet.add_observation(time, 'BJD_TDB', 30, 'mid', flux, np.full_like(flux, 0.001), 'flux', 'test')

    assert planet.observations[0]['observation_type'] == kind
    assert planet.observations[0]['epoch'] == 100


# ---------------------------------------------------------------------------
# 7: printing angles
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('angle, expected', [
    (plc.Degrees(10, 20, 30.05), '10:20:30.05'),
    (plc.Degrees(10, 20, 30.005), '10:20:30.005'),
    (plc.Degrees(10, 20, 30.5), '10:20:30.5'),
    (plc.Degrees(0, 0, 0.00001), '00:00:00.00001'),
    (plc.Degrees(10, 20, 59.9999999), '10:21:00.0'),
    (plc.Degrees(3, 15), '03:15:00.0'),
])
def test_07_dms_keeps_fractional_seconds(angle, expected):
    assert str(angle.dms()) == expected


def test_07_hms_and_coordinates_keep_fractional_seconds():
    assert str(plc.Hours(5, 6, 7.08).hms()) == '05:06:07.08'
    target = plc.FixedTarget(plc.Hours(5, 6, 7.08), plc.Degrees(-10, -20, -30.05))
    assert target.coord == '05:06:07.08 -10:20:30.05'


# ---------------------------------------------------------------------------
# 8: find_single_star on a non-square image
# ---------------------------------------------------------------------------

def _image_with_star(x, y, rows=60, columns=200):
    rng = np.random.default_rng(2)
    xx, yy = np.meshgrid(np.arange(columns) + 0.5, np.arange(rows) + 0.5)
    image = 100.0 + rng.normal(0, 5, (rows, columns))
    image += 2000.0 * np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / 2.0 ** 2)
    return image


def test_08_star_found_beyond_the_number_of_rows():
    image = _image_with_star(150.5, 30.5)

    star = plc.find_single_star(image, 150, 30, mean=100.0, std=5.0)

    assert star is not None
    assert star[0] == pytest.approx(150.5, abs=0.1)
    assert star[1] == pytest.approx(30.5, abs=0.1)


def test_08_predicted_y_outside_the_image():
    image = _image_with_star(150.5, 30.5)

    assert plc.find_single_star(image, 150, 80, mean=100.0, std=5.0) is None


# ---------------------------------------------------------------------------
# 10: exotethys model fallbacks (needs the databases)
# ---------------------------------------------------------------------------

def test_10_phoenix_irac_with_non_solar_metallicity_uses_atlas():
    ldcs = exotethys(4.5, 5500.0, 0.2, 'irac1', stellar_model='phoenix')
    assert len(ldcs) == 4 and np.all(np.isfinite(ldcs))


# ---------------------------------------------------------------------------
# 11: integer angles
# ---------------------------------------------------------------------------

def test_11_planet_accepts_integer_angles():
    as_int = plc.Planet('x', 330, 18, LOGG, TEFF, METALLICITY, RP, PERIOD, SMA, ECC, 90, 0, MID_TIME, 'BJD_TDB')
    as_float = plc.Planet('x', 330.0, 18.0, LOGG, TEFF, METALLICITY, RP, PERIOD, SMA, ECC, 90.0, 0.0, MID_TIME,
                          'BJD_TDB')

    assert as_int.inclination == as_float.inclination == 90.0
    assert as_int.periastron == as_float.periastron == 0.0
    assert as_int.target.ra.deg() == as_float.target.ra.deg()
    assert as_int.target.dec.deg() == as_float.target.dec.deg()


def test_11_locate_planet_accepts_integer_angles():
    planet = plc.locate_planet(331, 19, radius=plc.Degrees(2.0))
    assert isinstance(planet, plc.Planet)
