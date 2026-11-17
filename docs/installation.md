# Installation

PyLightcurve 4.0 is tested with Python 3.8 to 3.12.

## From PyPI

PyLightcurve 4.0 and the latest version are published under the same name,
`pylightcurve`. To install the latest 4.0 release, limit the version:

```bash
pip install "pylightcurve>=4.0,<4.1"
```

Plain `pip install pylightcurve` installs the latest version, PyLightcurve 5, which
is not compatible with code written for 4.0.

## From source

```bash
git clone https://github.com/atsiaras/pylightcurve-4.git
cd pylightcurve-4
pip install .
```

## Dependencies

pip installs these automatically:

| Package | Used for |
|:--|:--|
| [numpy](https://numpy.org), [scipy](https://scipy.org) | numerical calculations |
| [astropy](https://www.astropy.org) | FITS files |
| [emcee](https://emcee.readthedocs.io) | MCMC fitting |
| [matplotlib](https://matplotlib.org) | plots of the fitting results |
| [pyyaml](https://pyyaml.org) | catalogue files |
| [certifi](https://github.com/certifi/python-certifi) | secure downloads of the databases |

## Data downloaded on first use

The first time you `import pylightcurve`, it downloads the data it needs: the
limb-darkening grids, the filter passbands, the Exoplanet Characterisation
Catalogue and the ephemeris tables used for the time conversions. The data are
stored in `~/.pylightcurve4/`, and on later imports PyLightcurve checks them for
updates. The first import can therefore take a few minutes and needs an internet
connection.

If the data become corrupted, for example after an interrupted download, delete the
`~/.pylightcurve4/` folder and import PyLightcurve again.

## Check the installation

```python
import pylightcurve as plc

print(plc.__version__)
planet = plc.get_planet('HD209458b')
print(planet.transit_duration('COUSINS_R'))   # transit duration in days
```
