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
        self.repeated_magnitudes1 = [self.mag1, self.mag2]
        self.repeated_magnitudes3 = [self.mag4, self.mag5]

    def test_combine_magnitudes_equal_uncertainties(self):
        result = photometry.combine_observed_magnitudes(self.repeated_magnitudes1)
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

