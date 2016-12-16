import unittest
import numpy as np
from .context import superbol
import superbol.zip_photometry as zip_phot



class TestSortPhotometry(unittest.TestCase):
    
    def setUp(self):
        self.dtype = [('jd', '>f8'), ('name', 'S1'), ('magnitude', '>f8'), 
                      ('uncertainty', '>f8')]

    def test_sort_by_JD(self):
        """An unsorted list of photometry should be returned sorted by JD"""
        photometry = np.array([(2446855.622, 'U', 5.909, 0.02),
                               (2446853.583, 'B', 4.970, 0.02),
                               (2446854.635, 'U', 5.324, 0.02)], dtype=self.dtype)

        expected = np.array([(2446853.583, 'B', 4.970, 0.02),
                             (2446854.635, 'U', 5.324, 0.02),
                             (2446855.622, 'U', 5.909, 0.02)], dtype=self.dtype)
        result = zip_phot.sort_photometry(photometry)
        self.assertTrue(np.array_equal(expected, result))

class TestBinPhotometry(unittest.TestCase):
    
    def setUp(self):
        self.dtype = [('jd', '>f8'), ('name', 'S1'), ('magnitude', '>f8'), 
                      ('uncertainty', '>f8')]

    def test_bin_87A_March_1_data(self):
        """Data taken with one sunrise in between should land in two bins"""
        photometry = np.array([(2446855.622, 'U', 5.909, 0.02),
                               (2446855.622, 'B', 5.165, 0.02),
                               (2446855.622, 'V', 4.498, 0.02),
                               (2446855.622, 'R', 4.064, 0.02),
                               (2446855.622, 'I', 3.987, 0.02),
                               (2446856.270, 'J', 3.49, 0.02),
                               (2446856.270, 'H', 3.27, 0.02),
                               (2446856.270, 'K', 2.99, 0.02),
                               (2446856.270, 'L', 2.45, 0.03),
                               (2446856.270, 'M', 2.34, 0.04)], dtype=self.dtype)
        expected = np.array([(2446855.622, 'U', 5.909, 0.02),
                              (2446855.622, 'B', 5.165, 0.02),
                              (2446855.622, 'V', 4.498, 0.02),
                              (2446855.622, 'R', 4.064, 0.02),
                              (2446855.622, 'I', 3.987, 0.02),
                              (2446856.270, 'J', 3.49, 0.02),
                              (2446856.270, 'H', 3.27, 0.02),
                              (2446856.270, 'K', 2.99, 0.02),
                              (2446856.270, 'L', 2.45, 0.03),
                              (2446856.270, 'M', 2.34, 0.04)], dtype=self.dtype)
        result = zip_phot.bin_JDs(photometry)
        self.assertTrue(np.allclose(expected['jd'], result['jd']))

    def test_bin_87A_March_1_dt_1(self):
        """Using a large dt should group all the data from March 1"""
        photometry = np.array([(2446855.622, 'U', 5.909, 0.02),
                               (2446855.622, 'B', 5.165, 0.02),
                               (2446855.622, 'V', 4.498, 0.02),
                               (2446855.622, 'R', 4.064, 0.02),
                               (2446855.622, 'I', 3.987, 0.02),
                               (2446856.270, 'J', 3.49, 0.02),
                               (2446856.270, 'H', 3.27, 0.02),
                               (2446856.270, 'K', 2.99, 0.02),
                               (2446856.270, 'L', 2.45, 0.03),
                               (2446856.270, 'M', 2.34, 0.04)], dtype=self.dtype)
        
        expected = np.array([(2446855.946, 'U', 5.909, 0.02),
                             (2446855.946, 'B', 5.165, 0.02),
                             (2446855.946, 'V', 4.498, 0.02),
                             (2446855.946, 'R', 4.064, 0.02),
                             (2446855.946, 'I', 3.987, 0.02),
                             (2446855.946, 'J', 3.49, 0.02),
                             (2446855.946, 'H', 3.27, 0.02),
                             (2446855.946, 'K', 2.99, 0.02),
                             (2446855.946, 'L', 2.45, 0.03),
                             (2446855.946, 'M', 2.34, 0.04)], dtype=self.dtype)
        result = zip_phot.bin_JDs(photometry, dt=1.0)
        self.assertTrue(np.allclose(expected['jd'], result['jd']), msg='{0}, {1}'.format(expected, result))

    def test_bin_87A_March_2(self):
        """Using data from March 2 should force the merger of Vis and IR data to one date"""
        photometry = np.array([(2446856.639, 'U', 6.278, 0.02),
                               (2446856.639, 'B', 5.259, 0.02),
                               (2446856.639, 'V', 4.460, 0.02),
                               (2446856.639, 'R', 4.025, 0.02),
                               (2446856.639, 'I', 3.952, 0.02),
                               (2446856.510, 'J', 3.44, 0.02),
                               (2446856.510, 'H', 3.28, 0.02),
                               (2446856.510, 'K', 3.01, 0.02),
                               (2446856.510, 'L', 2.58, 0.03),
                               (2446856.510, 'M', 2.52, 0.04)], dtype=self.dtype)

        expected = np.array([(2446856.5745, 'U', 6.278, 0.02),
                             (2446856.5745, 'B', 5.259, 0.02),
                             (2446856.5745, 'V', 4.460, 0.02),
                             (2446856.5745, 'R', 4.025, 0.02),
                             (2446856.5745, 'I', 3.952, 0.02),
                             (2446856.5745, 'J', 3.44, 0.02),
                             (2446856.5745, 'H', 3.28, 0.02),
                             (2446856.5745, 'K', 3.01, 0.02),
                             (2446856.5745, 'L', 2.58, 0.03),
                             (2446856.5745, 'M', 2.52, 0.04)], dtype=self.dtype)

        result = zip_phot.bin_JDs(photometry)
        self.assertTrue(np.allclose(expected['jd'], result['jd']), msg='{0}, {1}'.format(expected, result))

class TestBinByDay(unittest.TestCase):
    
    def setUp(self):
        self.dtype = [('jd', '>f8'), ('name', 'S1'), ('magnitude', '>f8'), 
                      ('uncertainty', '>f8')]

    def test_bin_simple_dataset(self):
        photometry = np.array([(2446855.622, 'U', 5.909, 0.02),
                               (2446855.622, 'B', 5.165, 0.02),
                               (2446855.622, 'V', 4.498, 0.02),
                               (2446855.622, 'R', 4.064, 0.02),
                               (2446855.622, 'I', 3.987, 0.02),
                               (2446856.270, 'J', 3.49, 0.02),
                               (2446856.270, 'H', 3.27, 0.02),
                               (2446856.270, 'K', 2.99, 0.02),
                               (2446856.270, 'L', 2.45, 0.03),
                               (2446856.270, 'M', 2.34, 0.04)], dtype=self.dtype)
        expected = np.array([(2446855.0, 'U', 5.909, 0.02),
                              (2446855.0, 'B', 5.165, 0.02),
                              (2446855.0, 'V', 4.498, 0.02),
                              (2446855.0, 'R', 4.064, 0.02),
                              (2446855.0, 'I', 3.987, 0.02),
                              (2446856.0, 'J', 3.49, 0.02),
                              (2446856.0, 'H', 3.27, 0.02),
                              (2446856.0, 'K', 2.99, 0.02),
                              (2446856.0, 'L', 2.45, 0.03),
                              (2446856.0, 'M', 2.34, 0.04)], dtype=self.dtype)
        result = zip_phot.bin_by_day(photometry)
        self.assertTrue(np.allclose(expected['jd'], result['jd']))

class TestCombineRepeatedDetections(unittest.TestCase):

    def setUp(self):
        self.dtype = [('jd', '>f8'), ('name', 'S1'), ('magnitude', '>f8'), 
                      ('uncertainty', '>f8')]

    def test_combine_two_repeated_detections(self):
        """Two detections on the same JD with the same filter should be 
        averaged together"""
        err = np.sqrt(0.02**2 + 0.02**2)
        photometry = np.array([(2446856.01, 'J', 3.44, 0.02),
                               (2446856.01, 'J', 3.46, 0.02)], dtype=self.dtype)
        expected = np.array([(2446856.01, 'J', 3.45, err)], dtype=self.dtype)
        result = zip_phot.combine_repeated_detections(photometry)
        self.assertTrue(np.array_equal(expected, result), msg='{0}, {1}'.format(expected, result))

    def test_leave_non_repeated_detections_alone(self):
        """When there are no repeated detections, do nothing"""
        photometry = np.array([(2446856.01, 'J', 3.44, 0.02),
                               (2446856.01, 'H', 3.28, 0.02)], dtype=self.dtype)
        expected = photometry
        result = zip_phot.combine_repeated_detections(photometry)
        self.assertTrue(np.array_equal(expected, result), msg='{0}, {1}'.format(expected, result))

    def test_combine_repeated_detectiions_two_filters(self):
        photometry = np.array([(2446856.510, 'J', 3.44, 0.02),
                               (2446856.510, 'H', 3.28, 0.02),
                               (2446856.510, 'J', 3.46, 0.02),
                               (2446856.510, 'K', 3.01, 0.02), 
                               (2446856.510, 'L', 2.58, 0.03), 
                               (2446856.510, 'M', 2.52, 0.04),
                               (2446856.510, 'M', 2.55, 0.04),
                               (2446856.639, 'U', 6.278, 0.02), 
                               (2446856.639, 'B', 5.259, 0.02), 
                               (2446856.639, 'V', 4.460, 0.02), 
                               (2446856.639, 'R', 4.025, 0.02), 
                               (2446856.639, 'I', 3.952, 0.02)], dtype=self.dtype)

        jerr = np.sqrt(0.02**2 + 0.02**2)
        merr = np.sqrt(0.04**2 + 0.04**2)
        expected = np.array([(2446856.510, 'H', 3.28, 0.02),
                             (2446856.510, 'K', 3.01, 0.02), 
                             (2446856.510, 'L', 2.58, 0.03), 
                             (2446856.510, 'J', 3.45, jerr),
                             (2446856.510, 'M', 2.535, merr),
                             (2446856.639, 'U', 6.278, 0.02), 
                             (2446856.639, 'B', 5.259, 0.02), 
                             (2446856.639, 'V', 4.460, 0.02), 
                             (2446856.639, 'R', 4.025, 0.02), 
                             (2446856.639, 'I', 3.952, 0.02)], dtype=self.dtype)
        result = zip_phot.combine_repeated_detections(photometry)
        self.assertTrue(np.array_equal(expected, result), msg='{0}, {1}'.format(expected, result))
