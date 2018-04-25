import unittest
import numpy as np

from superbol import mag2flux
from superbol import photometry

class TestGroupObservedMagnitudes(unittest.TestCase):

    def setUp(self):

        self.mag10 = mag2flux.ObservedMagnitude(magnitude = 0,
                                                uncertainty = 0,
                                                band = 'U',
                                                time = 1)

        self.mag11 = mag2flux.ObservedMagnitude(magnitude = 0,
                                                uncertainty = 0,
                                                band = 'B',
                                                time = 1.1)

        self.mag12 = mag2flux.ObservedMagnitude(magnitude = 0,
                                                uncertainty = 0,
                                                band = 'V',
                                                time = 1.2)

        self.mag21 = mag2flux.ObservedMagnitude(magnitude = 0,
                                                uncertainty = 0,
                                                band = 'U',
                                                time = 2.1)

        self.mag27 = mag2flux.ObservedMagnitude(magnitude = 0,
                                                uncertainty = 0,
                                                band = 'B',
                                                time = 2.7)

    def test_group_magnitudes_floor_same_day(self):
        magnitudes = [self.mag10, self.mag11, self.mag12]
        expected = [[self.mag10, self.mag11, self.mag12]]
        result = photometry.group_magnitudes(magnitudes)
        self.assertEqual(expected, result)

    def test_group_magnitudes_floor_different_days(self):
        magnitudes = [self.mag10, self.mag11, self.mag12, self.mag21, self.mag27]
        expected = [[self.mag10, self.mag11, self.mag12], 
                    [self.mag21, self.mag27]]
        result = photometry.group_magnitudes(magnitudes)
        self.assertEqual(expected, result)

class TestCombineMagnitudes(unittest.TestCase):

    def setUp(self):
        self.mag1 = mag2flux.ObservedMagnitude(magnitude = 100,
                                               uncertainty = 10,
                                               band = 'U',
                                               time = 0)
        self.mag2 = mag2flux.ObservedMagnitude(magnitude = 200,
                                               uncertainty = 10,
                                               band = 'U',
                                               time = 0)
        self.mag3 = mag2flux.ObservedMagnitude(magnitude = 150,
                                               uncertainty = 8,
                                               band = 'B',
                                               time = 0)
        self.mag4 = mag2flux.ObservedMagnitude(magnitude = 50,
                                               uncertainty = 8,
                                               band = 'V',
                                               time = 0)
        self.mag5 = mag2flux.ObservedMagnitude(magnitude = 60,
                                               uncertainty = 8,
                                               band = 'V',
                                               time = 0)
        self.magnitudes = [self.mag1, self.mag2, self.mag3, self.mag4, self.mag5]
        self.repeated_magnitudes_U = [self.mag1, self.mag2]
        self.repeated_magnitudes_V = [self.mag4, self.mag5]

    def test_combine_magnitudes_equal_uncertainties(self):
        result = photometry.combine_observed_magnitudes(self.repeated_magnitudes_U)
        expected = mag2flux.ObservedMagnitude(magnitude = 150,
                                              uncertainty = np.sqrt(200)/2.,
                                              band = 'U',
                                              time = 0)
        self.assertEqual(expected, result)

    def test_combine_magnitudes_unequal_uncertainties(self):
        magnitudes = [self.mag1, self.mag3]
        result = photometry.combine_observed_magnitudes(magnitudes)
        expected = mag2flux.ObservedMagnitude(magnitude = 130.488,
                                              uncertainty = 6.247,
                                              band = 'U',
                                              time = 0)
        self.assertAlmostEqual(expected.magnitude, result.magnitude, 3)
        self.assertAlmostEqual(expected.uncertainty, result.uncertainty, 3)

    def test_yield_magnitudes_at_each_observed_band(self):
        result_generator = photometry.yield_observed_magnitudes_at_each_observed_band(self.magnitudes)
        expected_U = self.repeated_magnitudes_U
        expected_V = self.repeated_magnitudes_V
        self.assertEqual([self.mag3], next(result_generator))
        self.assertEqual(expected_U, next(result_generator))
        self.assertEqual(expected_V, next(result_generator))
        with self.assertRaises(StopIteration):
            next(result_generator)

    def test_get_multi_band_photometry(self):
        result = photometry.get_multi_band_photometry(self.magnitudes)
        expected = [self.mag3, mag2flux.ObservedMagnitude(150, np.sqrt(200)/2., 'U', 0),
                    mag2flux.ObservedMagnitude(55.0, np.sqrt(128)/2., 'V', 0)]

        # Ugly
        for i, flux in enumerate(result):
            self.assertAlmostEqual(result[i].magnitude, expected[i].magnitude)
            self.assertAlmostEqual(result[i].uncertainty, expected[i].uncertainty)

    def test_weighted_average(self):
        expected = 10.4
        result = photometry.weighted_average([10.0, 12.0], [0.5, 1.0])
        self.assertEqual(expected, result)

    def test_weighted_average_uncertainty(self):
        expected = 0.69
        uncertainties = [1, 1, 3]
        result = photometry.weighted_average_uncertainty(uncertainties)
        self.assertAlmostEqual(expected, result, 2)

    def test_get_weights(self):
        expected = [4, 1]
        uncertainties = [0.5, 1.0]
        result = photometry.get_weights(uncertainties)
        self.assertEqual(expected, result)

class TestInterpolateMultiBandPhotometry(unittest.TestCase):

    def setUp(self):
        self.mag01 = mag2flux.ObservedMagnitude(100, 2, 'U', 0)
        self.mag02 = mag2flux.ObservedMagnitude(200, 2, 'B', 0)
        self.mag03 = mag2flux.ObservedMagnitude(150, 2, 'V', 0)
        self.mag11 = mag2flux.ObservedMagnitude(100, 2, 'U', 1)
        self.mag13 = mag2flux.ObservedMagnitude(150, 2, 'V', 1)
        self.mag21 = mag2flux.ObservedMagnitude(100, 2, 'U', 2)
        self.mag22 = mag2flux.ObservedMagnitude(200, 2, 'B', 2)
        self.mag23 = mag2flux.ObservedMagnitude(150, 2, 'V', 2)
        self.mag31 = mag2flux.ObservedMagnitude(150, 2, 'U', 3)
        self.mag32 = mag2flux.ObservedMagnitude(150, 2, 'B', 3)
        self.multi_band_photometry0 = photometry.get_multi_band_photometry([self.mag01, self.mag02, self.mag03])
        self.multi_band_photometry1 = photometry.get_multi_band_photometry([self.mag11, self.mag13])
        self.multi_band_photometry2 = photometry.get_multi_band_photometry([self.mag21, self.mag22, self.mag23])
        self.multi_band_photometry3 = photometry.get_multi_band_photometry([self.mag31, self.mag32])

    def test_get_observed_band_names(self):
        expected = ['B', 'U', 'V']
        result = photometry.get_observed_band_names([self.multi_band_photometry0, self.multi_band_photometry1, self.multi_band_photometry2])
        self.assertEqual(expected, result)

    def test_get_observed_times(self):
        expected = [0, 1, 2]
        result = photometry.get_observed_times([self.multi_band_photometry0, self.multi_band_photometry1, self.multi_band_photometry2])
        self.assertEqual(expected, result)

    def test_get_unobserved_times(self):
        lightcurve = [self.mag02, self.mag22, self.mag32]
        observed_times = photometry.get_observed_times([self.multi_band_photometry0,
                                                        self.multi_band_photometry1,
                                                        self.multi_band_photometry2,
                                                        self.multi_band_photometry3])
        unobserved_times = photometry.get_unobserved_times(lightcurve, observed_times)
        self.assertEqual([1], unobserved_times)
