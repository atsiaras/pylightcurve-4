
__all__ = ['open_dict', 'save_dict', 'copy_dict',
           'open_yaml', 'save_yaml', 'copy_yaml',
           'open_fits', 'save_fits', 'copy_fits',
           'fits_sci', 'fits_err', 'fits_sci_err',
           'download'
           ]

import os
import ssl
import copy
import yaml
import numpy as np
import pickle
import urllib
import certifi

from urllib.request import urlretrieve
from astropy.io import fits as pf

def open_dict(path):
    """Load a dictionary from a pickle file.

    Files written by Python 2 are read with the ``latin-1`` encoding.

    Parameters
    ----------
    path : str
        Path of the pickle file.

    Returns
    -------
    object
        The unpickled object.

    """
    try:
        return pickle.load(open(path, 'rb'))
    except UnicodeDecodeError:
        return pickle.load(open(path, 'rb'), encoding='latin-1')


def save_dict(dictionary, path):
    """Save a dictionary to a pickle file (protocol 2).

    Parameters
    ----------
    dictionary : dict
        The dictionary to save.
    path : str
        Path of the pickle file.

    """
    internal_copy = copy_dict(dictionary)
    pickle.dump(dictionary, open(path, 'wb'), protocol=2)
    del internal_copy


def copy_dict(dictionary):
    """Return a deep copy of a dictionary.

    Parameters
    ----------
    dictionary : dict
        The dictionary to copy.

    Returns
    -------
    dict
        Independent copy.

    """
    return copy.deepcopy(dictionary)


def open_yaml(path):
    """Load a YAML file (safe loader).

    Parameters
    ----------
    path : str
        Path of the YAML file.

    Returns
    -------
    object
        The loaded content.

    """
    return yaml.load(open(path, 'r'), Loader=yaml.SafeLoader)


def save_yaml(dictionary, path):
    """Save a dictionary to a YAML file.

    Parameters
    ----------
    dictionary : dict
        The dictionary to save.
    path : str
        Path of the YAML file.

    """
    yaml.dump(dictionary, open(path, 'w'), default_flow_style=False)


def copy_yaml(dictionary):
    """Return a deep copy of a dictionary loaded from YAML.

    Parameters
    ----------
    dictionary : dict
        The dictionary to copy.

    Returns
    -------
    dict
        Independent copy.

    """
    return copy.deepcopy(dictionary)


def open_fits(path):
    """Open a FITS file and return an in-memory copy.

    The file is closed before returning, and header problems are fixed where
    possible.

    Parameters
    ----------
    path : str
        Path of the FITS file.

    Returns
    -------
    astropy.io.fits.HDUList
        Copy of the file content.

    """
    with pf.open(path, memmap=False) as hdulist:
        internal_copy = copy_fits(hdulist)
    internal_copy.verify('fix')
    return internal_copy


def save_fits(hdulist, path):
    """Save an HDU list to a FITS file, fixing header problems.

    Parameters
    ----------
    hdulist : astropy.io.fits.HDUList
        The data to save.
    path : str
        Path of the FITS file.

    """
    internal_copy = copy_fits(hdulist)
    try:
        internal_copy.writeto(path, output_verify='fix')
    except pf.VerifyError as e:
        hdu = int(str(e.args)[4:-4].split('\\n')[1].replace('HDU ', '').replace(':', ''))
        card = int(str(e.args)[4:-4].split('\\n')[2].replace('Card ', '').replace(':', ''))
        del internal_copy[hdu].header[card]
        internal_copy.writeto(path, output_verify='fix')
    
    del internal_copy


def copy_fits(hdulist):
    """Return a copy of an HDU list.

    Parameters
    ----------
    hdulist : astropy.io.fits.HDUList
        The HDU list to copy.

    Returns
    -------
    astropy.io.fits.HDUList
        Independent copy.

    """
    return pf.HDUList([ff.copy() for ff in hdulist])


def fits_sci(hdulist):
    """Find the indices of the ``SCI`` extensions of an HDU list.

    Parameters
    ----------
    hdulist : astropy.io.fits.HDUList
        The HDU list.

    Returns
    -------
    numpy.ndarray
        Indices of the extensions named ``SCI``.

    """
    return np.where(np.array([ff.name for ff in hdulist]) == 'SCI')[0]


def fits_err(hdulist):
    """Find the indices of the ``ERR`` extensions of an HDU list.

    Parameters
    ----------
    hdulist : astropy.io.fits.HDUList
        The HDU list.

    Returns
    -------
    numpy.ndarray
        Indices of the extensions named ``ERR``.

    """
    return np.where(np.array([ff.name for ff in hdulist]) == 'ERR')[0]


def fits_sci_err(hdulist):
    """Pair the ``SCI`` and ``ERR`` extensions of an HDU list.

    Parameters
    ----------
    hdulist : astropy.io.fits.HDUList
        The HDU list. It must have as many ``SCI`` as ``ERR`` extensions.

    Returns
    -------
    numpy.ndarray
        Array of ``[sci_index, err_index]`` pairs.

    """
    return np.swapaxes([fits_sci(hdulist), fits_err(hdulist)], 0, 1)


ctx = ssl.create_default_context(cafile=certifi.where())

def download(link, destination, filename=None):
    """Download a file.

    If the download fails, it is tried again with the certificate authorities
    provided by ``certifi``. The file is first written to ``destination + '.part'``
    and renamed when complete, so a failed download leaves no partial file.

    Parameters
    ----------
    link : str
        URL of the file.
    destination : str
        Local path to save it to.
    filename : str, optional
        Name shown in the progress messages. Default is the file name of
        ``destination``.

    Returns
    -------
    bool
        ``True`` if the download succeeded, ``False`` otherwise.

    """

    if not filename:
        filename = os.path.split(destination)[1]

    print('    Downloading {0}...'.format(filename))

    temporary = destination + '.part'

    try:
        urlretrieve(link, temporary)
        os.replace(temporary, destination)
        print('    Done!')
        return True
    except Exception:
        try:
            with urllib.request.urlopen(link, context=ctx) as u, open(temporary, 'wb') as f:
                f.write(u.read())
            os.replace(temporary, destination)
            print('    Done!')
            return True
        except Exception:
            try:
                os.remove(temporary)
            except OSError:
                pass
            print('    Could not download {0}'.format(filename))
            return False


def open_dict_online(link):

    try:
        return pickle.load(urllib.request.urlopen(link))
    except:
        try:
            return pickle.load(urllib.request.urlopen(link, context=ctx))
        except:
            return False
