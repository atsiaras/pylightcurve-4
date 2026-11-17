
import os
import pickle
import types

import numpy as np
import pytest

import pylightcurve.__databases__ as dbm
from pylightcurve.errors import PyLCLibraryError

ACCESSORS = {'exotethys': 'exotethys', 'ephemerides': 'ephemeris',
             'photometry': 'photometry', 'catalogues': 'catalogues'}

DATABASES = ['exotethys', 'ephemerides', 'photometry', 'catalogues']

VERSION = '.'.join(open(os.path.join(os.path.dirname(dbm.__file__), '__version__.txt')).read().split('.')[:2])


class FakeRemote:
    """A local folder that plays the part of the online servers."""

    def __init__(self, root, home):
        self.root = root
        self.home = home
        self.downloads = []      # every link requested, in order
        self.index = {}          # database name -> link of its current pickle

    def put(self, link, content):
        """Publish a file. ``content`` is bytes, str, or any picklable object."""
        path = os.path.join(self.root, link)
        if isinstance(content, str):
            content = content.encode()
        elif not isinstance(content, bytes):
            content = pickle.dumps(content)
        with open(path, 'wb') as f:
            f.write(content)

    def publish(self, name, files, frequency=0, tag='v1', zipfile_link='no_zip.zip'):
        """Publish a database.

        ``files`` maps a key to ``(local_path, link)`` or
        ``(local_path, link, external_link)``. The file contents must be put()
        separately (or deliberately left missing).
        """
        entries = {}
        for key, spec in files.items():
            entries[key] = {'local_path': spec[0], 'link': spec[1]}
            if len(spec) > 2:
                entries[key]['external_link'] = spec[2]
        link = '{0}_{1}.pickle'.format(name, tag)
        self.put(link, {'files': entries, 'frequency': frequency, 'zipfile': zipfile_link})
        self.index[name] = link
        self._write_index()

    def publish_all_empty(self, skip=()):
        for name in DATABASES:
            if name not in skip and name not in self.index:
                self.publish(name, {})

    def _write_index(self):
        directory = os.path.join(self.home, '.' + dbm.package_name)
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, dbm.databases_file), 'wb') as f:
            pickle.dump({VERSION: dict(self.index)}, f)

    def local(self, *parts):
        return os.path.join(self.home, '.' + dbm.package_name, *parts)

    def download(self, link, destination, filename=None):
        self.downloads.append(link)
        source = os.path.join(self.root, link)
        if not os.path.isfile(source):
            return False
        with open(source, 'rb') as fin, open(destination, 'wb') as fout:
            fout.write(fin.read())
        return True


@pytest.fixture
def remote(tmp_path, monkeypatch):
    home = tmp_path / 'home'
    root = tmp_path / 'remote'
    home.mkdir()
    root.mkdir()
    monkeypatch.setenv('HOME', str(home))
    fake = FakeRemote(str(root), str(home))
    monkeypatch.setattr(dbm, 'download', fake.download)
    monkeypatch.setattr(dbm, 'open_dict_online', lambda link: False)
    return fake


def freeze_today(monkeypatch, date, hour='12:00:00'):
    """Pretend that today is ``date`` ('YYYY-MM-DD') for the update-frequency check.

    Freezes both clocks the code reads: ``astrotime.now()`` (deciding if an update
    is due) and ``time.strftime`` (writing the date of the last update).
    """
    now = '{0} {1}'.format(date, hour)
    monkeypatch.setattr(dbm.astrotime, 'now', classmethod(lambda cls: cls(now)))
    yymmdd = date[2:4] + date[5:7] + date[8:10]
    monkeypatch.setattr(dbm, 'time', types.SimpleNamespace(strftime=lambda fmt: yymmdd))



def test_first_install_downloads_everything(remote):
    remote.put('stars.yaml', 'a: 1\n')
    remote.put('planets.yaml', 'b: 2\n')
    remote.publish('catalogues', {'stars': ('catalogues/ecc_stars.yaml', 'stars.yaml'),
                                  'planets': ('catalogues/ecc_planets.yaml', 'planets.yaml')})
    remote.publish_all_empty()

    data = dbm.PlcData()

    for name in DATABASES:
        assert getattr(data, name if name != 'ephemerides' else 'ephemeris')() == remote.local(name)
    assert open(remote.local('catalogues', 'ecc_stars.yaml')).read() == 'a: 1\n'
    assert open(remote.local('catalogues_link.txt')).read() == 'catalogues_v1.pickle'
    assert data.ecc() == {'stars': {'a': 1}, 'planets': {'b': 2}}


def test_second_start_downloads_nothing(remote):
    remote.put('stars.yaml', 'a: 1\n')
    remote.publish('catalogues', {'stars': ('catalogues/ecc_stars.yaml', 'stars.yaml')})
    remote.publish_all_empty()
    dbm.PlcData()

    remote.downloads.clear()
    dbm.PlcData()

    assert remote.downloads == []


def test_changed_link(remote):
    remote.put('stars_v1.yaml', 'version: 1\n')
    remote.put('planets.yaml', 'b: 2\n')
    remote.publish('catalogues', {'stars': ('catalogues/ecc_stars.yaml', 'stars_v1.yaml'),
                                  'planets': ('catalogues/ecc_planets.yaml', 'planets.yaml')})
    remote.publish_all_empty()
    dbm.PlcData()

    remote.put('stars_v2.yaml', 'version: 2\n')
    remote.publish('catalogues', {'stars': ('catalogues/ecc_stars.yaml', 'stars_v2.yaml'),
                                  'planets': ('catalogues/ecc_planets.yaml', 'planets.yaml')}, tag='v2')
    remote.downloads.clear()
    data = dbm.PlcData()

    assert open(remote.local('catalogues', 'ecc_stars.yaml')).read() == 'version: 2\n'
    assert open(remote.local('catalogues_link.txt')).read() == 'catalogues_v2.pickle'
    assert 'planets.yaml' not in remote.downloads
    assert data.ecc()['stars'] == {'version': 2}


def test_no_update(remote):
    remote.put('stars.yaml', 'a: 1\n')
    remote.publish('catalogues', {'stars': ('catalogues/ecc_stars.yaml', 'stars.yaml')})
    remote.publish_all_empty()
    dbm.PlcData()

    remote.index['catalogues'] = 'catalogues_v2_not_on_server.pickle'
    remote._write_index()
    data = dbm.PlcData()

    assert data.catalogues() == remote.local('catalogues')
    assert open(remote.local('catalogues', 'ecc_stars.yaml')).read() == 'a: 1\n'


def test_deleted_local_file_is_downloaded_again(remote):
    remote.put('stars.yaml', 'a: 1\n')
    remote.publish('catalogues', {'stars': ('catalogues/ecc_stars.yaml', 'stars.yaml')})
    remote.publish_all_empty()
    dbm.PlcData()

    os.remove(remote.local('catalogues', 'ecc_stars.yaml'))
    dbm.PlcData()

    assert open(remote.local('catalogues', 'ecc_stars.yaml')).read() == 'a: 1\n'


def test_all_filters_and_ecc_are_cached(remote, monkeypatch):
    remote.put('stars.yaml', 'a: 1\n')
    remote.put('planets.yaml', 'b: 2\n')
    remote.put('filter.pickle', {'x': 1})
    remote.publish('catalogues', {'stars': ('catalogues/ecc_stars.yaml', 'stars.yaml'),
                                  'planets': ('catalogues/ecc_planets.yaml', 'planets.yaml')})
    remote.publish('photometry', {'r': ('photometry/COUSINS_R.pickle', 'filter.pickle'),
                                  'v': ('photometry/JOHNSON_V.pickle', 'filter.pickle')})
    remote.publish_all_empty()
    data = dbm.PlcData()

    yaml_calls = []
    real_open_yaml = dbm.open_yaml
    monkeypatch.setattr(dbm, 'open_yaml', lambda path: yaml_calls.append(path) or real_open_yaml(path))

    assert sorted(data.all_filters()) == ['COUSINS_R', 'JOHNSON_V']
    data.ecc()
    data.ecc()
    assert len(yaml_calls) == 2

    os.remove(remote.local('photometry', 'JOHNSON_V.pickle'))
    assert sorted(data.all_filters()) == ['COUSINS_R', 'JOHNSON_V']


def _ephemeris_table(jd):
    t = jd + np.linspace(-0.5, 1.5, 21)
    return {jd: {'t': t, 'ra': 2.0 * (t - jd) ** 2 + 1.0, 'dec': 3.0 * (t - jd) - 5.0,
                 'd': 1.0 + 0.01 * (t - jd), 'dt': 0.5 * (t - jd) ** 3}}


@pytest.mark.parametrize('method, filename', [('barycentre', 'bjd_dict.pickle'),
                                              ('heliocentre', 'hjd_dict.pickle')])
def test_ephemeris_interpolation(remote, method, filename):
    jd = 2460000
    remote.put('table.pickle', _ephemeris_table(jd))
    remote.publish('ephemerides', {'t': ('ephemerides/' + filename, 'table.pickle')})
    remote.publish_all_empty()
    data = dbm.PlcData()

    x = 0.37
    ra, dec, d, dt = getattr(data, method)(jd + x)

    # cubic interpolation is exact for polynomials up to third degree
    assert ra == pytest.approx(2.0 * x ** 2 + 1.0, rel=1e-9)
    assert dec == pytest.approx(3.0 * x - 5.0, rel=1e-9)
    assert d == pytest.approx(1.0 + 0.01 * x, rel=1e-9)
    assert dt == pytest.approx(0.5 * x ** 3, rel=1e-9)

    with pytest.raises(KeyError):
        getattr(data, method)(jd + 5.2)


def test_file_dropped_from_database_is_deleted(remote):
    remote.put('a.txt', 'a')
    remote.put('b.txt', 'b')
    remote.publish('exotethys', {'a': ('exotethys/a.txt', 'a.txt'), 'b': ('exotethys/b.txt', 'b.txt')})
    remote.publish_all_empty()
    dbm.PlcData()
    assert os.path.isfile(remote.local('exotethys', 'b.txt'))

    remote.publish('exotethys', {'a': ('exotethys/a.txt', 'a.txt')}, tag='v2')
    dbm.PlcData()

    assert os.path.isfile(remote.local('exotethys', 'a.txt'))
    assert not os.path.isfile(remote.local('exotethys', 'b.txt'))


def test_failed_unzip(remote, monkeypatch):
    class BrokenZip:
        def __init__(self, *args, **kwargs):
            pass

        def extractall(self, *args, **kwargs):
            raise OSError('corrupt archive')

    monkeypatch.setattr(dbm, 'zipfile', types.SimpleNamespace(ZipFile=BrokenZip))
    remote.put('archive.zip', b'not really a zip')
    remote.publish('exotethys', {}, zipfile_link='archive.zip')
    remote.publish_all_empty()
    monkeypatch.chdir(remote.root)

    dbm.PlcData()

    assert os.getcwd() == remote.root


def test_unavailable_database(remote):
    remote.publish('catalogues', {'stars': ('catalogues/ecc_stars.yaml', 'not_on_server.yaml')})
    remote.publish_all_empty()
    data = dbm.PlcData()

    with pytest.raises(PyLCLibraryError):
        data.ecc()


def test_unknown_version(remote):
    remote.publish_all_empty()
    directory = remote.local()
    with open(os.path.join(directory, dbm.databases_file), 'wb') as f:
        pickle.dump({'0.0': dict(remote.index)}, f)

    with pytest.raises(PyLCLibraryError):
        dbm.PlcData()


@pytest.mark.parametrize('broken', DATABASES)
def test_missing_file_disables_only_that_database(remote, capsys, broken):
    remote.publish(broken, {'f': (broken + '/f.txt', 'not_on_server.txt')})
    remote.publish_all_empty()

    data = dbm.PlcData()

    assert '{0} features cannot be used'.format(broken) in capsys.readouterr().out
    with pytest.raises(PyLCLibraryError, match='not available'):
        getattr(data, ACCESSORS[broken])()
    for name in DATABASES:
        if name != broken:
            assert getattr(data, ACCESSORS[name])() == remote.local(name)


def test_missing_database_pickle_disables_that_database(remote):
    remote.publish_all_empty(skip=['photometry'])
    remote.index['photometry'] = 'photometry_not_on_server.pickle'
    remote._write_index()

    data = dbm.PlcData()

    with pytest.raises(PyLCLibraryError):
        data.photometry()
    with pytest.raises(PyLCLibraryError):
        data.all_filters()
    assert data.catalogues() == remote.local('catalogues')


def test_database_recovers_when_download_works_again(remote):
    remote.publish('catalogues', {'stars': ('catalogues/ecc_stars.yaml', 'stars.yaml')})
    remote.publish_all_empty()
    data = dbm.PlcData()
    with pytest.raises(PyLCLibraryError):
        data.catalogues()

    remote.put('stars.yaml', 'a: 1\n')
    data = dbm.PlcData()

    assert data.catalogues() == remote.local('catalogues')


def test_builtin_index_copied_on_first_run(remote):
    local_index = remote.local(dbm.databases_file)
    assert not os.path.exists(local_index)

    dbm.PlcData()

    with open(local_index, 'rb') as a, open(os.path.join(os.path.dirname(dbm.__file__), dbm.databases_file), 'rb') as b:
        assert a.read() == b.read()


def test_newer_online_index_is_saved_and_used(remote, monkeypatch):
    remote.put('stars_v1.yaml', 'version: 1\n')
    remote.publish('catalogues', {'stars': ('catalogues/ecc_stars.yaml', 'stars_v1.yaml')})
    remote.publish_all_empty()
    dbm.PlcData()

    # the server publishes v2, but only the online index knows about it
    local_index = open(remote.local(dbm.databases_file), 'rb').read()
    remote.put('stars_v2.yaml', 'version: 2\n')
    remote.publish('catalogues', {'stars': ('catalogues/ecc_stars.yaml', 'stars_v2.yaml')}, tag='v2')
    online_index = {VERSION: dict(remote.index)}
    with open(remote.local(dbm.databases_file), 'wb') as f:
        f.write(local_index)
    monkeypatch.setattr(dbm, 'open_dict_online', lambda link: online_index)

    data = dbm.PlcData()

    assert data.databases == online_index
    assert pickle.load(open(remote.local(dbm.databases_file), 'rb')) == online_index
    assert open(remote.local('catalogues', 'ecc_stars.yaml')).read() == 'version: 2\n'



def publish_with_external(remote, frequency=7):
    remote.put('eop.txt', 'old')
    remote.put('eop_external.txt', 'fresh')
    remote.publish('ephemerides', {'eop': ('ephemerides/eop.txt', 'eop.txt', 'eop_external.txt')},
                   frequency=frequency)
    remote.publish_all_empty()


def test_first_install_updates_external_links(remote, monkeypatch):
    publish_with_external(remote)
    freeze_today(monkeypatch, '2026-11-10')

    data = dbm.PlcData()

    assert 'eop_external.txt' in remote.downloads
    assert open(remote.local('ephemerides', 'eop.txt')).read() == 'fresh'
    assert open(remote.local('ephemerides_last_update.txt')).read() == '261110'
    assert data.ephemeris() == remote.local('ephemerides')


@pytest.mark.parametrize('last_update, today, expected', [
    ('2026-11-10', '2026-11-10', False),   # same day
    ('2026-11-10', '2026-11-16', False),   # 6 days later
    ('2026-11-10', '2026-11-17', True),    # 7 days later
    ('2026-11-28', '2026-12-01', False),   # 3 days later, across a month
    ('2026-11-28', '2026-12-05', True),    # 7 days later, across a month
    ('2026-12-28', '2027-01-01', False),   # 4 days later, across a year
    ('2026-12-28', '2027-01-04', True),    # 7 days later, across a year
    ('2028-02-25', '2028-03-01', False),   # 5 days later, across a leap day
    ('2028-02-25', '2028-03-03', True),    # 7 days later, across a leap day
])
def test_external_links_follow_update_frequency(remote, monkeypatch, last_update, today, expected):
    publish_with_external(remote)
    freeze_today(monkeypatch, last_update)
    dbm.PlcData()

    remote.downloads.clear()
    freeze_today(monkeypatch, today)
    dbm.PlcData()

    assert ('eop_external.txt' in remote.downloads) is expected
    expected_date = today if expected else last_update
    assert open(remote.local('ephemerides_last_update.txt')).read() == expected_date[2:4] + expected_date[5:7] + expected_date[8:10]


def test_frequency_zero_never_updates_external_links(remote, monkeypatch):
    publish_with_external(remote, frequency=0)
    freeze_today(monkeypatch, '2026-11-10')

    dbm.PlcData()

    assert 'eop_external.txt' not in remote.downloads
    assert open(remote.local('ephemerides', 'eop.txt')).read() == 'old'
    assert not os.path.exists(remote.local('ephemerides_last_update.txt'))


def test_unreadable_last_update_file_triggers_update(remote, monkeypatch):
    publish_with_external(remote)
    freeze_today(monkeypatch, '2026-11-10')
    dbm.PlcData()
    with open(remote.local('ephemerides_last_update.txt'), 'w') as f:
        f.write('')

    remote.downloads.clear()
    freeze_today(monkeypatch, '2026-11-11')
    dbm.PlcData()

    assert 'eop_external.txt' in remote.downloads


def test_missing_local_file_resets_last_update(remote, monkeypatch):
    publish_with_external(remote)
    freeze_today(monkeypatch, '2026-11-10')
    dbm.PlcData()
    os.remove(remote.local('ephemerides', 'eop.txt'))

    remote.downloads.clear()
    freeze_today(monkeypatch, '2026-11-11')
    dbm.PlcData()

    # the file is downloaded again and, since the last-update record was reset,
    # refreshed from its external link straight away
    assert remote.downloads.index('eop.txt') < remote.downloads.index('eop_external.txt')
    assert open(remote.local('ephemerides', 'eop.txt')).read() == 'fresh'

def test_failed_external_update_keeps_local_copy_usable(remote, monkeypatch):
    publish_with_external(remote)
    freeze_today(monkeypatch, '2026-11-10')
    dbm.PlcData()

    os.remove(os.path.join(remote.root, 'eop_external.txt'))   # external server down
    freeze_today(monkeypatch, '2026-11-20')
    data = dbm.PlcData()

    assert data.ephemeris() == remote.local('ephemerides')
    assert open(remote.local('ephemerides', 'eop.txt')).read() == 'fresh'
