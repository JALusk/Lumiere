import unittest

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

