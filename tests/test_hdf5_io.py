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

