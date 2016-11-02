import unittest
import tables as tb
from .context import superbol
import superbol.data.hdf5_io as hdf5_io

class TestNewSNGroupCreation(unittest.TestCase):

    def setUp(self):
        self.hdf5_filename = 'tests/test_data.h5'
        self.h5file = tb.open_file(self.hdf5_filename, 'a')
        self.sn_name = 'test_sn'

    def test_make_new_sn_group_if_sn_group_exists(self):
        result = hdf5_io.make_new_sn_group(self.h5file, 'sn1998a')
        expected = None
        self.assertEqual(expected, result)

    def test_make_new_sn_group_if_sn_group_does_not_exist(self):
        hdf5_io.make_new_sn_group(self.h5file, self.sn_name)
        group_exists = self.h5file.__contains__('/sn/' + self.sn_name)
        self.assertTrue(group_exists)
        if group_exists:
            self.h5file.remove_node('/sn', self.sn_name, recursive = True)

    def tearDown(self):
        self.h5file.close()

class TestNewPhotTableCreation(unittest.TestCase):

    def setUp(self):
        self.hdf5_filename = 'tests/test_data.h5'
        self.h5file = tb.open_file(self.hdf5_filename, 'a')
        self.sn_name = 'test_sn'

    def test_make_new_sn_phot_table_if_table_exists(self):
        result = hdf5_io.make_new_sn_phot_table(self.h5file, 'sn1998a')
        expected = self.h5file.get_node('/sn/sn1998a/phot')
        self.assertEqual(expected, result)

    def test_make_new_sn_phot_table_if_table_does_not_exist(self):
        hdf5_io.make_new_sn_group(self.h5file, self.sn_name)
        phot_table = hdf5_io.make_new_sn_phot_table(self.h5file, self.sn_name)
        table_exists = self.h5file.__contains__('/sn/' + self.sn_name + '/phot')
        self.assertTrue(table_exists)
        if table_exists:
            self.h5file.remove_node('/sn', self.sn_name, recursive = True)

    def tearDown(self):
        self.h5file.close()

class TestNewParametersTableCreation(unittest.TestCase):

    def setUp(self):
        self.hdf5_filename = 'tests/test_data.h5'
        self.h5file = tb.open_file(self.hdf5_filename, 'a')
        self.sn_name = 'test_sn'

    def test_make_new_sn_parameters_table_if_table_exists(self):
        result = hdf5_io.make_new_sn_parameters_table(self.h5file, 'sn1998a')
        expected = self.h5file.get_node('/sn/sn1998a/parameters')
        self.assertEqual(expected, result)

    def test_make_new_sn_parameters_table_if_table_does_not_exist(self):
        hdf5_io.make_new_sn_group(self.h5file, self.sn_name)
        phot_table = hdf5_io.make_new_sn_parameters_table(self.h5file, self.sn_name)
        table_exists = self.h5file.__contains__('/sn/' + self.sn_name + '/parameters')
        self.assertTrue(table_exists)
        if table_exists:
            self.h5file.remove_node('/sn', self.sn_name, recursive = True)

    def tearDown(self):
        self.h5file.close()

class TestNewFilter(unittest.TestCase):

    def setUp(self):
        self.hdf5_filename = 'tests/test_data.h5'
        self.h5file = tb.open_file(self.hdf5_filename, 'a')
        self.filter_table = self.h5file.root.filters
        self.filter_name = 'U'
        self.filter_eff_wl = 3660.0
        self.filter_flux_zeropoint = 4.175E-9
        self.note = 'TEST'
        self.ref = 'TEST'

    def test_set_new_filter_id(self):
        result = hdf5_io.set_new_filter_id(self.filter_table)
        expected = 99
        self.assertEqual(expected, result)

    def test_get_old_filter_id(self):
        result = hdf5_io.get_old_filter_id(self.filter_table, self.filter_name)
        expected = 4
        self.assertEqual(expected, result)

    def tearDown(self):
        self.h5file.close()

class TestMakeNewObservationEntry(unittest.TestCase):

    def setUp(self):
        self.hdf5_filename = 'tests/test_data.h5'
        self.h5file = tb.open_file(self.hdf5_filename, 'a')
        self.phot_table = self.h5file.root.sn.sn1998a.phot
        self.filter_id = 1
        self.jd = 9999999
        self.mag = 17.0
        self.err = 1.0
        self.ref = 'TESTREF'
        self.note = 'TESTNOTE'

    def test_make_new_observation_returns_row(self):
        new_observation = hdf5_io.make_new_observation_entry(self.filter_id, self.jd, self.mag, self.err, self.ref, self.note, self.phot_table)
        result = isinstance(new_observation, tb.tableextension.Row)
        self.assertTrue(result)

    def test_make_new_observation_returns_row_with_correct_filter_id(self):
        new_observation = hdf5_io.make_new_observation_entry(self.filter_id, self.jd, self.mag, self.err, self.ref, self.note, self.phot_table)
        expected = self.filter_id
        result = new_observation['filter_id']
        self.assertEqual(expected, result)

    def test_make_new_observation_returns_row_with_correct_jd(self):
        new_observation = hdf5_io.make_new_observation_entry(self.filter_id, self.jd, self.mag, self.err, self.ref, self.note, self.phot_table)
        expected = self.jd
        result = new_observation['jd']
        self.assertEqual(expected, result)

    def test_make_new_observation_returns_row_with_correct_mag(self):
        new_observation = hdf5_io.make_new_observation_entry(self.filter_id, self.jd, self.mag, self.err, self.ref, self.note, self.phot_table)
        expected = self.mag
        result = new_observation['magnitude']
        self.assertEqual(expected, result)

    def test_make_new_observation_returns_row_with_correct_err(self):
        new_observation = hdf5_io.make_new_observation_entry(self.filter_id, self.jd, self.mag, self.err, self.ref, self.note, self.phot_table)
        expected = self.err
        result = new_observation['uncertainty']
        self.assertEqual(expected, result)

    def test_make_new_observation_returns_row_with_correct_ref(self):
        new_observation = hdf5_io.make_new_observation_entry(self.filter_id, self.jd, self.mag, self.err, self.ref, self.note, self.phot_table)
        expected = self.ref
        result = new_observation['reference']
        self.assertEqual(expected, result)

    def test_make_new_observation_returns_row_with_correct_note(self):
        new_observation = hdf5_io.make_new_observation_entry(self.filter_id, self.jd, self.mag, self.err, self.ref, self.note, self.phot_table)
        expected = self.note
        result = new_observation['note']
        self.assertEqual(expected, result)

    def test_make_new_observation_returns_row_with_expected_telescope_id(self):
        new_observation = hdf5_io.make_new_observation_entry(self.filter_id, self.jd, self.mag, self.err, self.ref, self.note, self.phot_table)
        expected = 0
        result = new_observation['telescope_id']
        self.assertEqual(expected, result)

    def tearDown(self):
        self.h5file.close()
