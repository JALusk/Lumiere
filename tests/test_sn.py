import unittest
from .context import superbol
import superbol.sn as sn
import tables as tb
import numpy as np
from pkg_resources import resource_filename

class TestSNInitialization(unittest.TestCase):
   
    def setUp(self):
        self.sn_name = "SN1987A"
        self.my_sn = sn.SN(self.sn_name)

    def test_sn_name_initialization(self):
        expected = self.sn_name
        result = self.my_sn.name
        self.assertEqual(expected, result)

    def test_sn_source_initialization_when_no_source_given(self):
        expected = None
        result = self.my_sn.source
        self.assertEqual(expected, result)

    def test_sn_source_initialization_when_source_given(self):
        source = '/path/to/source'
        my_sn = sn.SN(self.sn_name, source)
        expected = source
        result = my_sn.source
        self.assertEqual(expected, result)

    def test_sn_filter_table_initialization(self):
        expected = None
        result = self.my_sn.filter_table
        self.assertEqual(expected, result)

    def test_sn_phot_table_initialization(self):
        expected = None
        result = self.my_sn.phot_table
        self.assertEqual(expected, result)

    def test_sn_parameter_table_initialization(self):
        expected = None
        result = self.my_sn.parameter_table
        self.assertEqual(expected, result)

class TestSNDefaultSourceOpening(unittest.TestCase):

    def setUp(self):
        self.sn_name = "sn1998a"
        self.default_source = None
        self.my_sn = sn.SN(self.sn_name, self.default_source)

    def test_sn_open_default_source_file_returns_pytables_file_object(self):
        h5file = self.my_sn.open_source_h5file()
        result = isinstance(h5file, tb.file.File)
        h5file.close()
        self.assertTrue(result)

    def test_sn_open_default_source_file_returns_default_sn_data_file(self):
        h5file = self.my_sn.open_source_h5file()
        expected = resource_filename('superbol', 'data/sn_data.h5')
        result = h5file.filename
        h5file.close()
        self.assertEqual(expected, result)

class TestSNCustomSourceOpening(unittest.TestCase):

    def setUp(self):
        self.sn_name = "sn1998a"
        self.custom_source = "tests/test_data.h5"
        self.my_sn = sn.SN(self.sn_name, self.custom_source)

    def test_sn_open_custom_source_file_returns_pytables_file_object(self):
        h5file = self.my_sn.open_source_h5file()
        result = isinstance(h5file, tb.file.File)
        h5file.close()
        self.assertTrue(result)

    def test_sn_open_custom_source_file_returns_correct_filename(self):
        h5file = self.my_sn.open_source_h5file()
        expected = self.custom_source
        result = h5file.filename
        h5file.close()
        self.assertEqual(expected, result)

class TestSNImportHDF5Tables(unittest.TestCase):

    def setUp(self):
        self.sn_name = "sn1998a"
        self.source = "tests/test_data.h5"
        self.my_sn = sn.SN(self.sn_name, self.source)
        self.h5file = tb.open_file(self.source, 'r')
    
    def test_import_hdf5_tables_sets_filter_table_to_group_object(self):
        self.my_sn.import_hdf5_tables(self.h5file)
        result = isinstance(self.my_sn.filter_table, tb.Table)
        self.assertTrue(result)

    def test_import_hdf5_tables_sets_correct_filter_table(self):
        self.my_sn.import_hdf5_tables(self.h5file)
        expected = self.h5file.root.filters
        result = self.my_sn.filter_table
        self.assertEqual(expected, result)

    def test_import_hdf5_tables_sets_phot_table_to_group_object(self):
        self.my_sn.import_hdf5_tables(self.h5file)
        result = isinstance(self.my_sn.phot_table, tb.Table)
        self.assertTrue(result)

    def test_import_hdf5_tables_sets_correct_phot_table(self):
        self.my_sn.import_hdf5_tables(self.h5file)
        expected = self.h5file.get_node("/sn/"+self.sn_name, "phot")
        result = self.my_sn.phot_table
        self.assertEqual(expected, result)

    def test_import_hdf5_tables_sets_parameter_table_to_group_object(self):
        self.my_sn.import_hdf5_tables(self.h5file)
        result = isinstance(self.my_sn.parameter_table, tb.Table)
        self.assertTrue(result)

    def test_import_hdf5_tables_sets_correct_parameter_table(self):
        self.my_sn.import_hdf5_tables(self.h5file)
        expected = self.h5file.get_node("/sn/"+self.sn_name, "parameters")
        result = self.my_sn.parameter_table
        self.assertEqual(expected, result)

    def tearDown(self):
        self.h5file.close()

class TestSNGetPhotometry(unittest.TestCase):

    def setUp(self):
        self.sn_name = "sn1998a"
        self.source = "tests/test_data.h5"
        self.my_sn = sn.SN(self.sn_name, self.source)
        self.h5file = tb.open_file(self.source, 'r')
        self.my_sn.import_hdf5_tables(self.h5file)
    
    def test_get_photometry_returns_numpy_array(self):
        photometry = self.my_sn.get_photometry()
        result = isinstance(photometry, np.ndarray)
        self.assertTrue(result)

    def test_get_photometry_dtype(self):
        expected = [('jd', '>f8'), ('name', 'S1'), ('magnitude', '>f8'), ('uncertainty', '>f8')]
        photometry = self.my_sn.get_photometry()
        result = photometry.dtype
        self.assertEqual(expected, result)

    def test_get_photometry_98A_first_obs_jd(self):
        expected = 2450820.3
        
        photometry = self.my_sn.get_photometry()
        result = photometry[0]['jd']
        self.assertEqual(expected, result)

    def test_get_photometry_98A_first_obs_filter_name(self):
        expected = b'R'
        
        photometry = self.my_sn.get_photometry()
        result = photometry[0]['name']
        self.assertEqual(expected, result)

    def test_get_photometry_98A_first_obs_mag(self):
        expected = 17.0
        
        photometry = self.my_sn.get_photometry()
        result = photometry[0]['magnitude']
        self.assertEqual(expected, result)

    def test_get_photometry_98A_first_obs_uncertainty(self):
        expected = 0.5
        
        photometry = self.my_sn.get_photometry()
        result = photometry[0]['uncertainty']
        self.assertEqual(expected, result)

    def tearDown(self):
        self.h5file.close()

class TestSNGetBCColor(unittest.TestCase):

    def setUp(self):
        self.sn_name = "sn1998a"
        self.source = "tests/test_data.h5"
        self.my_sn = sn.SN(self.sn_name, self.source)
        self.h5file = tb.open_file(self.source, 'r')
        self.my_sn.import_hdf5_tables(self.h5file)
        self.photometry = self.my_sn.get_photometry()

    def test_get_color_returns_float_when_given_valid_data(self):
        color = self.my_sn.get_color(self.photometry, 2450837.8, 'B', 'V')
        result = isinstance(color, float)
        self.assertTrue(result)

    def test_get_color_returns_None_when_given_bad_JD(self):
        color = self.my_sn.get_color(self.photometry, 12, 'B', 'V')
        self.assertEqual(color, None)

    def test_get_color_returns_None_when_given_bad_filter1(self):
        color = self.my_sn.get_color(self.photometry, 2450837.8, 'X', 'V')
        self.assertEqual(color, None)

    def test_get_color_returns_None_when_given_bad_filter2(self):
        color = self.my_sn.get_color(self.photometry, 2450837.8, 'B', 'X')
        self.assertEqual(color, None)

    def test_get_color_98A_returns_first_BV_color_correctly(self):
        B = 18.05
        V = 16.92
        expected = B-V
        result = self.my_sn.get_color(self.photometry, 2450837.8, 'B', 'V')
        self.assertEqual(expected, result)
 
    def tearDown(self):
        self.h5file.close()

class TestSNGetBCColorUncertainty(unittest.TestCase):

    def setUp(self):
        self.sn_name = "sn1998a"
        self.source = "tests/test_data.h5"
        self.my_sn = sn.SN(self.sn_name, self.source)
        self.h5file = tb.open_file(self.source, 'r')
        self.my_sn.import_hdf5_tables(self.h5file)
        self.photometry = self.my_sn.get_photometry()

    def test_get_color_uncertainty_returns_float_when_given_valid_data(self):
        color = self.my_sn.get_color_uncertainty(self.photometry, 2450837.8, 'B', 'V')
        result = isinstance(color, float)
        self.assertTrue(result)

    def test_get_color_returns_None_when_given_bad_JD(self):
        color = self.my_sn.get_color_uncertainty(self.photometry, 12, 'B', 'V')
        self.assertEqual(color, None)

    def test_get_color_returns_None_when_given_bad_filter1(self):
        color = self.my_sn.get_color_uncertainty(self.photometry, 2450837.8, 'X', 'V')
        self.assertEqual(color, None)

    def test_get_color_returns_None_when_given_bad_filter2(self):
        color = self.my_sn.get_color_uncertainty(self.photometry, 2450837.8, 'B', 'X')
        self.assertEqual(color, None)

    def test_get_color_98A_returns_first_BV_color_correctly(self):
        B_err = 0.1
        V_err = 0.03
        expected = np.sqrt(B_err**2 + V_err**2)
        result = self.my_sn.get_color_uncertainty(self.photometry, 2450837.8, 'B', 'V')
        self.assertEqual(expected, result)
 
    def tearDown(self):
        self.h5file.close()

class TestSNDereddenUBVRIMagnitudes(unittest.TestCase):

    def setUp(self):
        self.sn_name = "sn1998a"
        self.source = "tests/test_data.h5"
        self.my_sn = sn.SN(self.sn_name, self.source)
        self.h5file = tb.open_file(self.source, 'r')
        self.my_sn.import_hdf5_tables(self.h5file)
        self.photometry = self.my_sn.get_photometry()

    def test_deredden_UBVRI_magnitudes_98A_works_on_first_R_magnitude(self):
        Av_gal = self.my_sn.parameter_table.cols.Av_gal[0]
        Av_host = self.my_sn.parameter_table.cols.Av_host[0]
        Av_tot = Av_gal + Av_host
        expected = 17.0 - 0.751 * Av_tot
        
        dereddened_photometry = self.my_sn.deredden_UBVRI_magnitudes(self.photometry)
        result = dereddened_photometry[0]['magnitude']

        self.assertEqual(expected, result)

    def test_deredden_UBVRI_magnitudes_98A_ignores_J_magnitude(self):
        expected = 15.17
        
        dereddened_photometry = self.my_sn.deredden_UBVRI_magnitudes(self.photometry)
        result = dereddened_photometry[20]['magnitude']

        self.assertEqual(expected, result)

    def tearDown(self):
        self.h5file.close()
