"""More tests for pylightcurve/processes/files.py.

These cover the parts that test_processes_files.py does not reach or does not check:
the Latin-1 fallback of open_dict, every path of download (including the retry with
certifi and the clean-up of partial downloads), open_dict_online, and the results of
the FITS helpers. Everything runs offline: downloads use local file:// URLs.
"""

import os
import pickle
import ssl

import numpy as np
import pytest
from astropy.io import fits as pf

import pylightcurve as plc
import pylightcurve.processes.files as files


def file_url(path):
    return 'file://' + os.path.abspath(str(path))


@pytest.fixture
def source(tmp_path):
    path = tmp_path / 'source.txt'
    path.write_text('new content')
    return path


# ---------------------------------------------------------------------------
# open_dict / save_dict
# ---------------------------------------------------------------------------

def test_open_dict_reads_python2_pickles(tmp_path):
    # written by Python 2: pickle.dumps({'name': 'h\xe9'}, protocol=2)
    path = tmp_path / 'python2.pickle'
    path.write_bytes(b'\x80\x02}q\x00U\x04nameq\x01U\x02h\xe9q\x02s.')

    assert plc.open_dict(str(path)) == {'name': 'h\xe9'}


def test_save_dict_round_trip_does_not_change_the_input(tmp_path):
    data = {'a': [1, 2, 3], 'b': {'c': np.arange(3)}}
    path = str(tmp_path / 'data.pickle')

    plc.save_dict(data, path)
    loaded = plc.open_dict(path)

    assert loaded['a'] == [1, 2, 3]
    assert np.array_equal(loaded['b']['c'], np.arange(3))
    assert data['a'] == [1, 2, 3]


def test_copy_dict_and_copy_yaml_are_deep_copies():
    for copy_function in [plc.copy_dict, plc.copy_yaml]:
        data = {'a': [1, 2]}
        copied = copy_function(data)
        copied['a'].append(3)
        assert data == {'a': [1, 2]}


# ---------------------------------------------------------------------------
# download
# ---------------------------------------------------------------------------

def test_download_success(tmp_path, source, capsys):
    destination = tmp_path / 'downloaded.txt'

    assert plc.download(file_url(source), str(destination)) is True

    assert destination.read_text() == 'new content'
    assert not os.path.exists(str(destination) + '.part')
    output = capsys.readouterr().out
    assert 'Downloading downloaded.txt' in output and 'Done!' in output


def test_download_uses_the_given_name_in_messages(tmp_path, source, capsys):
    plc.download(file_url(source), str(tmp_path / 'downloaded.txt'), filename='my database')

    assert 'Downloading my database' in capsys.readouterr().out


def test_download_failure_creates_nothing(tmp_path, capsys):
    destination = tmp_path / 'downloaded.txt'

    assert plc.download(file_url(tmp_path / 'missing.txt'), str(destination)) is False

    assert not destination.exists()
    assert not os.path.exists(str(destination) + '.part')
    assert 'Could not download downloaded.txt' in capsys.readouterr().out


def test_download_failure_keeps_existing_file(tmp_path):
    destination = tmp_path / 'downloaded.txt'
    destination.write_text('old content')

    assert plc.download(file_url(tmp_path / 'missing.txt'), str(destination)) is False

    assert destination.read_text() == 'old content'


def test_download_retries_with_certifi_when_first_attempt_fails(tmp_path, source, monkeypatch):
    def failing_urlretrieve(link, destination):
        raise OSError('certificate verify failed')

    monkeypatch.setattr(files, 'urlretrieve', failing_urlretrieve)
    destination = tmp_path / 'downloaded.txt'

    assert plc.download(file_url(source), str(destination)) is True

    assert destination.read_text() == 'new content'
    assert not os.path.exists(str(destination) + '.part')


def test_interrupted_download_leaves_no_partial_file(tmp_path, monkeypatch):
    def interrupted_urlretrieve(link, destination):
        with open(destination, 'w') as f:
            f.write('half of the fi')
        raise OSError('connection reset')

    monkeypatch.setattr(files, 'urlretrieve', interrupted_urlretrieve)
    destination = tmp_path / 'downloaded.txt'
    destination.write_text('old content')

    assert plc.download(file_url(tmp_path / 'missing.txt'), str(destination)) is False

    assert destination.read_text() == 'old content'
    assert not os.path.exists(str(destination) + '.part')


def test_retry_context_verifies_certificates():
    assert files.ctx.verify_mode == ssl.CERT_REQUIRED
    assert files.ctx.check_hostname is True


# ---------------------------------------------------------------------------
# open_dict_online
# ---------------------------------------------------------------------------

def test_open_dict_online_success(tmp_path):
    path = tmp_path / 'online.pickle'
    path.write_bytes(pickle.dumps({'version': 2}))

    assert files.open_dict_online(file_url(path)) == {'version': 2}


def test_open_dict_online_failure_returns_false(tmp_path):
    assert files.open_dict_online(file_url(tmp_path / 'missing.pickle')) is False


def test_open_dict_online_retries_with_certifi(tmp_path, monkeypatch):
    path = tmp_path / 'online.pickle'
    path.write_bytes(pickle.dumps({'version': 2}))

    real_urlopen = files.urllib.request.urlopen
    calls = []

    def urlopen(link, context=None):
        calls.append(context)
        if context is None:
            raise OSError('certificate verify failed')
        return real_urlopen(link, context=context)

    monkeypatch.setattr(files.urllib.request, 'urlopen', urlopen)

    assert files.open_dict_online(file_url(path)) == {'version': 2}
    assert calls == [None, files.ctx]


# ---------------------------------------------------------------------------
# FITS helpers (need astropy)
# ---------------------------------------------------------------------------

def make_hdulist():
    return pf.HDUList([pf.PrimaryHDU(),
                       pf.ImageHDU(data=np.array([[1.0, 2.0], [3.0, 4.0]]), name='SCI'),
                       pf.ImageHDU(data=np.array([[0.1, 0.2], [0.3, 0.4]]), name='ERR'),
                       pf.ImageHDU(data=np.array([[5.0, 6.0], [7.0, 8.0]]), name='SCI'),
                       pf.ImageHDU(data=np.array([[0.5, 0.6], [0.7, 0.8]]), name='ERR')])


def test_fits_sci_err_with_several_extensions():
    hdulist = make_hdulist()

    assert list(plc.fits_sci(hdulist)) == [1, 3]
    assert list(plc.fits_err(hdulist)) == [2, 4]
    assert plc.fits_sci_err(hdulist).tolist() == [[1, 2], [3, 4]]


def test_fits_sci_err_without_extensions():
    hdulist = pf.HDUList([pf.PrimaryHDU()])

    assert len(plc.fits_sci(hdulist)) == 0
    assert len(plc.fits_err(hdulist)) == 0


def test_copy_fits_is_independent():
    original = make_hdulist()
    copied = plc.copy_fits(original)

    copied[1].data[0, 0] = 99.0
    copied[1].header['NEW'] = 1

    assert original[1].data[0, 0] == 1.0
    assert 'NEW' not in original[1].header


def test_save_and_open_fits_round_trip(tmp_path):
    path = str(tmp_path / 'image.fits')
    original = make_hdulist()
    original[1].header['EXPTIME'] = 30.0

    plc.save_fits(original, path)
    loaded = plc.open_fits(path)
    os.remove(path)  # open_fits returns an in-memory copy

    assert [hdu.name for hdu in loaded] == ['PRIMARY', 'SCI', 'ERR', 'SCI', 'ERR']
    assert np.array_equal(loaded[3].data, original[3].data)
    assert loaded[1].header['EXPTIME'] == 30.0


def add_invalid_card(header):
    # Card.fromstring does not validate, so the illegal keyword 'A?X' reaches the header;
    # astropy refuses to write it, which sends save_fits into its repair branch
    header.append(pf.Card.fromstring('A?X     =                    3'.ljust(80)))


def has_invalid_card(header):
    return any(key.upper() == 'A?X' for key in header.keys())


def test_save_fits_removes_invalid_card_and_keeps_the_rest(tmp_path):
    path = str(tmp_path / 'image.fits')
    original = make_hdulist()
    original[1].header['EXPTIME'] = 30.0
    add_invalid_card(original[1].header)
    assert has_invalid_card(original[1].header)

    plc.save_fits(original, path)
    loaded = plc.open_fits(path)

    assert not has_invalid_card(loaded[1].header)
    assert loaded[1].header['EXPTIME'] == 30.0
    assert np.array_equal(loaded[1].data, original[1].data)
    # the input is not changed
    assert has_invalid_card(original[1].header)
