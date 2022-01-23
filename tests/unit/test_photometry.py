import math
import unittest
import numpy as np

from superbol import mag2flux
from superbol import photometry

U_band = mag2flux.Band('U', 'CTIO U', 3660.0, 417.5E-11)
B_band = mag2flux.Band('B', 'CTIO B', 4380.0, 632.0E-11)
V_band = mag2flux.Band('V', 'CTIO V', 5450.0, 363.1E-11)

class TestGroupObservedMagnitudes(unittest.TestCase):

    def setUp(self):

        self.mag10 = mag2flux.ObservedMagnitude(magnitude = 0,
                                                uncertainty = 0,
                                                band = U_band,
                                                time = 1)

        self.mag11 = mag2flux.ObservedMagnitude(magnitude = 0,
                                                uncertainty = 0,
                                                band = B_band,
                                                time = 1.1)

        self.mag12 = mag2flux.ObservedMagnitude(magnitude = 0,
                                                uncertainty = 0,
                                                band = V_band,
                                                time = 1.2)

        self.mag21 = mag2flux.ObservedMagnitude(magnitude = 0,
                                                uncertainty = 0,
                                                band = U_band,
                                                time = 2.1)

        self.mag27 = mag2flux.ObservedMagnitude(magnitude = 0,
                                                uncertainty = 0,
                                                band = B_band,
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
                                               band = U_band,
                                               time = 0)
        self.mag2 = mag2flux.ObservedMagnitude(magnitude = 200,
                                               uncertainty = 10,
                                               band = U_band,
                                               time = 0)
        self.mag3 = mag2flux.ObservedMagnitude(magnitude = 150,
                                               uncertainty = 8,
                                               band = B_band,
                                               time = 0)
        self.mag4 = mag2flux.ObservedMagnitude(magnitude = 50,
                                               uncertainty = 8,
                                               band = V_band,
                                               time = 0)
        self.mag5 = mag2flux.ObservedMagnitude(magnitude = 60,
                                               uncertainty = 8,
                                               band = V_band,
                                               time = 0)
        self.magnitudes = [self.mag1, self.mag2, self.mag3, self.mag4, self.mag5]
        self.repeated_magnitudes_U = [self.mag1, self.mag2]
        self.repeated_magnitudes_V = [self.mag4, self.mag5]

    def test_combine_magnitudes_equal_uncertainties(self):
        result = photometry.combine_observed_magnitudes(self.repeated_magnitudes_U)
        expected = mag2flux.ObservedMagnitude(magnitude = 150,
                                              uncertainty = np.sqrt(200)/2.,
                                              band = U_band,
                                              time = 0)
        self.assertEqual(expected, result)

    def test_combine_magnitudes_unequal_uncertainties(self):
        magnitudes = [self.mag1, self.mag3]
        result = photometry.combine_observed_magnitudes(magnitudes)
        expected = mag2flux.ObservedMagnitude(magnitude = 130.488,
                                              uncertainty = 6.247,
                                              band = U_band,
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
        expected = [self.mag3, mag2flux.ObservedMagnitude(150, np.sqrt(200)/2., U_band, 0),
                    mag2flux.ObservedMagnitude(55.0, np.sqrt(128)/2., V_band, 0)]

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
        self.mag01 = mag2flux.ObservedMagnitude(100, 2, U_band, 0)
        self.mag02 = mag2flux.ObservedMagnitude(200, 2, B_band, 0)
        self.mag03 = mag2flux.ObservedMagnitude(150, 2, V_band, 0)
        self.mag11 = mag2flux.ObservedMagnitude(100, 2, U_band, 1)
        self.mag13 = mag2flux.ObservedMagnitude(150, 2, V_band, 1)
        self.mag21 = mag2flux.ObservedMagnitude(100, 2, U_band, 2)
        self.mag22 = mag2flux.ObservedMagnitude(200, 2, B_band, 2)
        self.mag23 = mag2flux.ObservedMagnitude(150, 2, V_band, 2)
        self.mag31 = mag2flux.ObservedMagnitude(150, 2, U_band, 3)
        self.mag32 = mag2flux.ObservedMagnitude(150, 2, B_band, 3)
        self.multi_band_photometry0 = photometry.get_multi_band_photometry([self.mag01, self.mag02, self.mag03])
        self.multi_band_photometry1 = photometry.get_multi_band_photometry([self.mag11, self.mag13])
        self.multi_band_photometry2 = photometry.get_multi_band_photometry([self.mag21, self.mag22, self.mag23])
        self.multi_band_photometry3 = photometry.get_multi_band_photometry([self.mag31, self.mag32])

    # TODO check
    def test_get_observed_band_names(self):
        expected = ['B', 'U', 'V']
        result = photometry.get_observed_band_names([self.multi_band_photometry0, self.multi_band_photometry1, self.multi_band_photometry2])
        self.assertEqual(expected, result)

    def test_get_times(self):
        expected = [0, 1, 2, 3]
        result = photometry.get_times([self.mag01, self.mag13, self.mag21, self.mag32])
        self.assertEqual(expected, result)

    # TODO check
    def test_get_observed_times(self):
        expected = [0, 1, 2]
        result = photometry.get_observed_times([self.multi_band_photometry0, self.multi_band_photometry1, self.multi_band_photometry2])
        self.assertEqual(expected, result)

    # TODO check
    def test_get_unobserved_times(self):
        lightcurve = [self.mag02, self.mag22, self.mag32]
        observed_times = photometry.get_observed_times([self.multi_band_photometry0,
                                                        self.multi_band_photometry1,
                                                        self.multi_band_photometry2,
                                                        self.multi_band_photometry3])
        unobserved_times = photometry.get_unobserved_times(lightcurve, observed_times)
        self.assertEqual([1], unobserved_times)

    #check
    def test_get_lightcurve(self):
        expected = [self.mag01, self.mag11, self.mag21]
        result = photometry.get_lightcurve([self.multi_band_photometry0,
                                            self.multi_band_photometry1,
                                            self.multi_band_photometry2], 'U')
        self.assertEqual(expected, result)

    #check
    def test_get_previous_observed_magnitude(self):
        lightcurve = [self.mag02, self.mag22, self.mag32]
        unobserved_time = 1
        previous_observed_magnitude = photometry.get_previous_observed_magnitude(lightcurve, unobserved_time)
        self.assertEqual(self.mag02, previous_observed_magnitude)

    def test_get_previous_observed_magnitude_out_of_bounds(self):
        lightcurve = [self.mag22, self.mag32]
        unobserved_time = 1
        with self.assertRaises(photometry.MissingMagnitudeOutOfBounds):
            photometry.get_previous_observed_magnitude(lightcurve, unobserved_time)

    #check
    def test_get_next_observed_magnitude(self):
        lightcurve = [self.mag02, self.mag22, self.mag32]
        unobserved_time = 1
        next_observed_magnitude = photometry.get_next_observed_magnitude(lightcurve, unobserved_time)
        self.assertEqual(self.mag22, next_observed_magnitude)

    def test_get_next_observed_magnitude_out_of_bounds(self):
        lightcurve = [self.mag02, self.mag22, self.mag32]
        unobserved_time = 4
        with self.assertRaises(photometry.MissingMagnitudeOutOfBounds):
            photometry.get_next_observed_magnitude(lightcurve, unobserved_time)

    #check
    def test_get_interpolated_magnitudes(self):
        lightcurve = [self.mag02, self.mag22, self.mag32]
        previous_observed_magnitude = self.mag02
        next_observed_magnitude = self.mag22
        unobserved_time = 1
        weight1 = (2-1)/(2-0)
        weight2 = (1-0)/(2-0)
        uncertainty = math.sqrt(weight1**2 * previous_observed_magnitude.uncertainty**2 + weight2**2 + next_observed_magnitude.uncertainty**2)
        expected = [mag2flux.ObservedMagnitude(200.0, uncertainty , B_band, 1)]
        result = photometry.get_interpolated_magnitudes(lightcurve, [1])
        self.assertEqual(expected, result)

    #check
    def test_get_interpolated_magnitude_uncertainty(self):
        previous_magnitude = self.mag02
        next_magnitude = self.mag22
        unobserved_time = 1
        weight1 = (2-1)/(2-0)
        weight2 = (1-0)/(2-0)
        expected = math.sqrt(weight1**2 * previous_magnitude.uncertainty**2 + weight2**2 + next_magnitude.uncertainty**2)
        result = photometry.get_interpolated_magnitude_uncertainty(previous_magnitude, next_magnitude, unobserved_time)
        self.assertEqual(expected, result)


    #check
    def test_big_gap(self):
        lightcurve = [self.mag01, self.mag31]
        interpolated_lc = photometry.get_interpolated_magnitudes(lightcurve, [0, 1, 2, 3])
        self.assertEqual([], interpolated_lc)


    #check
    def test_append_interpolated_magnitudes_to_multi_band_photometry_set(self):
        magnitude = mag2flux.ObservedMagnitude(200.0, 2, 'B', 1)
        multi_band_photometry = self.multi_band_photometry1
        photometry.append_interpolated_magnitudes_to_multi_band_photometry_set([magnitude], [multi_band_photometry])
        self.assertTrue(magnitude in multi_band_photometry)

    #check
    def test_simple_interpolation(self):
        photometry.interpolate_missing_magnitudes([self.multi_band_photometry0,
                                                   self.multi_band_photometry1,
                                                   self.multi_band_photometry2])
        previous_magnitude = self.mag02
        next_magnitude = self.mag22
        unobserved_time = 1
        weight1 = (2-1)/(2-0)
        weight2 = (1-0)/(2-0)
        uncertainty = math.sqrt(weight1**2 * previous_magnitude.uncertainty**2 + weight2**2 + next_magnitude.uncertainty**2)
        interpolated_magnitude = mag2flux.ObservedMagnitude(200, uncertainty, B_band, 1)
        self.assertTrue(interpolated_magnitude in self.multi_band_photometry1)
