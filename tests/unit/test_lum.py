import unittest
import math

from superbol import lum

class TestDistance(unittest.TestCase):

    def setUp(self):
        self.distance = lum.Distance(1, 0.1)

    def test_distance_value(self):
        expected = 1
        result = self.distance.value
        self.assertEqual(expected, result)

    def test_distance_uncertainty(self):
        expected = 0.1
        result = self.distance.uncertainty
        self.assertEqual(expected, result)

    def test_negative_distance_raises_exception(self):
        with self.assertRaises(ValueError):
            negative_distance = lum.Distance(-1, 0.1)

class TestBolometricLuminosity(unittest.TestCase):

    def setUp(self):
        self.lbol = lum.BolometricLuminosity(1E42, 1E41)

    def test_lum_value(self):
        expected = 1E42
        result = self.lbol.value
        self.assertEqual(expected, result)

    def test_distance_uncertainty(self):
        expected = 1E41
        result = self.lbol.uncertainty
        self.assertEqual(expected, result)

    def test_negative_distance_raises_exception(self):
        with self.assertRaises(ValueError):
            negative_lbol = lum.BolometricLuminosity(-1, 0.1)

class TestBolometricFlux(unittest.TestCase):

    def setUp(self):
        self.fbol = lum.BolometricFlux(2, 0.2)

    def test_lum_value(self):
        expected = 2
        result = self.fbol.value
        self.assertEqual(expected, result)

    def test_distance_uncertainty(self):
        expected = 0.2
        result = self.fbol.uncertainty
        self.assertEqual(expected, result)

    def test_negative_distance_raises_exception(self):
        with self.assertRaises(ValueError):
            negative_fbol = lum.BolometricFlux(-1, 0.1)

    def test_convert_flux_value_to_lum_value(self):
        distance = lum.Distance(2, 0.2)
        expected = 4.0 * math.pi * 4 * 2
        result = self.fbol.to_lbol(distance)
        self.assertEqual(expected, result.value)

    def test_convert_flux_uncertainty_to_lum_uncertainty(self):
        distance = lum.Distance(2, 0.2)
        expected = 22.479
        result = self.fbol.to_lbol(distance)
        self.assertAlmostEqual(expected, result.uncertainty, 3)

class TestConvertFluxToLuminosity(unittest.TestCase):

    def test_convert_zero_flux(self):
        expected = 0
        fbol = lum.BolometricFlux(0, 0)
        distance = lum.Distance(1, 0.1)
        result = lum.convert_flux_to_luminosity(fbol, distance)
        self.assertEqual(expected, result.value)

    def test_convert_flux_value_to_lum_value(self):
        expected = 22.479
        fbol = lum.BolometricFlux(2, 0.2)
        distance = lum.Distance(2, 0.2)
        result = lum.convert_flux_to_luminosity(fbol, distance)
        self.assertAlmostEqual(expected, result.uncertainty, 3)

