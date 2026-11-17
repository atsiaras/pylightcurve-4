import pytest
import numpy as np
import pylightcurve as plc
import os
from astropy.io import fits as pf


def add_invalid_card(header):
    # Card.fromstring does not validate, so the illegal keyword 'A?X' reaches the header;
    # astropy refuses to write it, which sends save_fits into its repair branch
    header.append(pf.Card.fromstring('A?X     =                    3'.ljust(80)))


def has_invalid_card(header):
    return any(key.upper() == 'A?X' for key in header.keys())


def test_files():

    __location__ = os.path.abspath(os.path.dirname(__file__))

    xx = {'a': 1}
    plc.save_dict(xx, os.path.join(__location__, 'test.pickle'))
    xx = plc.open_dict(os.path.join(__location__, 'test.pickle'))
    os.remove(os.path.join(__location__, 'test.pickle'))
    yy = plc.copy_dict(xx)
    assert id(yy) != id(xx)
    del xx, yy

    xx = {'a': 1}
    plc.save_yaml(xx, os.path.join(__location__, 'test.yaml'))
    xx = plc.open_yaml(os.path.join(__location__, 'test.yaml'))
    os.remove(os.path.join(__location__, 'test.yaml'))
    yy = plc.copy_yaml(xx)
    assert id(yy) != id(xx)
    del xx, yy

    if os.path.isfile(os.path.join(__location__, 'test.fits')):
        os.remove(os.path.join(__location__, 'test.fits'))
    xx = pf.HDUList([pf.PrimaryHDU(), pf.ImageHDU(data=np.array([[1, 2], [1, 2]]), name='SCI'),
                     pf.ImageHDU(data=np.array([[1, 2], [1, 2]]), name='ERR')])

    print(plc.fits_sci(xx), plc.fits_err(xx), plc.fits_sci_err(xx))

    assert plc.fits_sci(xx)[0] == 1
    assert plc.fits_err(xx)[0] == 2
    assert plc.fits_sci_err(xx)[0][0] == 1
    assert plc.fits_sci_err(xx)[0][1] == 2

    add_invalid_card(xx[1].header)
    assert has_invalid_card(xx[1].header)

    plc.save_fits(xx, os.path.join(__location__, 'test.fits'))
    xx = plc.open_fits(os.path.join(__location__, 'test.fits'))
    os.remove(os.path.join(__location__, 'test.fits'))
    yy = plc.copy_fits(xx)
    assert id(yy) != id(xx)
    del xx, yy

    plc.download('https://raw.githubusercontent.com/ucl-exoplanets/pylightcurve/master/README.md',
                 os.path.join(__location__, 'test.txt'))
    os.remove(os.path.join(__location__, 'test.txt'))

    plc.download('https://raw.githubusercontent.com/ucl-exoplanets/pylightcurve/master/README.mx',
                 os.path.join(__location__, 'test.txt'))


