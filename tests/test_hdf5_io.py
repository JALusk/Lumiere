import unittest
import tables as tb
import numpy as np
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
        self.jd = 9999999.00
        self.mag = 17.0
        self.err = 1.0
        self.ref = b'TESTREF'
        self.note = b'TESTNOTE'

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

class TestAddSNFilterPhotometry(unittest.TestCase):

    def setUp(self):
        self.hdf5_filename = 'tests/test_data.h5'
        self.sn_name = 'test_sn'

    def test_add_sn_filter_photometry_old_filter(self):
        filter_name = 'U'
        filter_id = 4
        filter_eff_wl = None
        filter_flux_zeropoint = None
        ref = 'TESTREF'
        note = 'TESTNOTE'
        jd = [9999999.00, 10000000.00]
        mag = [17.0, 18.0]
        err = [0.1, 0.2]

        hdf5_io.add_sn_filter_photometry(self.hdf5_filename, self.sn_name, filter_name, filter_eff_wl, filter_flux_zeropoint, ref, note, jd, mag, err)

        h5file = tb.open_file(self.hdf5_filename, 'a')
        phot_table = h5file.root.sn.test_sn.phot
    
        dtype = [('filter_id', '<i8'), ('jd', '<f8'), ('magnitude', '<f8'), ('note', 'S32'), ('reference', 'S32'), ('telescope_id', '<i8'), ('uncertainty', '<f8')]
        expected = np.array([(filter_id, jd[0], mag[0], note, ref, 0, err[0])], dtype=dtype)
        result = phot_table[0]
        np.testing.assert_array_equal(expected, result)
        
        h5file.remove_node('/sn', self.sn_name, recursive = True)

        h5file.close()

    def test_add_sn_filter_photometry_new_filter(self):
        filter_name = 'Q'
        filter_id = 99
        filter_eff_wl = 3000.0
        filter_flux_zeropoint = 4.175E-9
        ref = 'TESTREF'
        note = 'TESTNOTE'
        jd = [9999999.00, 10000000.00]
        mag = [17.0, 18.0]
        err = [0.1, 0.2]

        hdf5_io.add_sn_filter_photometry(self.hdf5_filename, self.sn_name, filter_name, filter_eff_wl, filter_flux_zeropoint, ref, note, jd, mag, err)

        h5file = tb.open_file(self.hdf5_filename, 'a')
        phot_table = h5file.root.sn.test_sn.phot
    
        dtype = [('filter_id', '<i8'), ('jd', '<f8'), ('magnitude', '<f8'), ('note', 'S32'), ('reference', 'S32'), ('telescope_id', '<i8'), ('uncertainty', '<f8')]
        expected = np.array([(filter_id, jd[0], mag[0], note, ref, 0, err[0])], dtype=dtype)
        result = phot_table[0]
        np.testing.assert_array_equal(expected, result)
        
        h5file.remove_node('/sn', self.sn_name, recursive = True)
        h5file.root.filters.remove_rows(-1)

        h5file.close()


class TestAddSNParameters(unittest.TestCase):

    def setUp(self):
        self.hdf5_filename = 'tests/test_data.h5'
        self.sn_name = 'test_sn'

    def test_add_sn_parameters(self):
        Av_gal = 0.399
        Av_gal_ref = '1998ApJ...500..525S'
        Av_host = 0.0
        Av_host_ref = '2005MNRAS.360..950P'
        distance_Mpc = 30.34
        distance_Mpc_err = 7.0
        distance_Mpc_ref = '2005MNRAS.360..950P'
        explosion_JD = 2450801.0
        explosion_JD_err = 4.0
        explosion_JD_ref = '2005MNRAS.360..950P'
        heliocentric_v_kms = 2090.0
        heliocentric_v_kms_err = 2.0
        heliocentric_v_kms_ref = '2004AJ....128...16K'


        hdf5_io.add_sn_parameters(self.hdf5_filename, self.sn_name, Av_gal, Av_gal_ref, Av_host, Av_host_ref, distance_Mpc, distance_Mpc_err, distance_Mpc_ref, explosion_JD, explosion_JD_err, explosion_JD_ref, heliocentric_v_kms, heliocentric_v_kms_err, heliocentric_v_kms_ref)

        h5file = tb.open_file(self.hdf5_filename, 'a')
        parameters_table = h5file.root.sn.test_sn.parameters

        dtype = [('Av_gal', '<f8'), ('Av_gal_ref', 'S19'), ('Av_host', '<f8'), ('Av_host_ref', 'S19'), ('distance_Mpc', '<f8'), ('distance_Mpc_err', '<f8'), ('distance_Mpc_ref', 'S19'), ('explosion_JD', '<f8'), ('explosion_JD_err', '<f8'), ('explosion_JD_ref', 'S19'), ('heliocentric_v_kms', '<f8'), ('heliocentric_v_kms_err', '<f8'), ('heliocentric_v_kms_ref', 'S19')]
        expected = np.array([(Av_gal, Av_gal_ref, Av_host, Av_host_ref, distance_Mpc, distance_Mpc_err, distance_Mpc_ref, explosion_JD, explosion_JD_err, explosion_JD_ref, heliocentric_v_kms, heliocentric_v_kms_err, heliocentric_v_kms_ref)], dtype=dtype)
        result = parameters_table[0]
        np.testing.assert_array_equal(expected, result)

        h5file.remove_node('/sn', self.sn_name, recursive = True)

        h5file.close()

