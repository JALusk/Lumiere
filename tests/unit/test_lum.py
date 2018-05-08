import unittest
import math

from superbol import luminosity

class TestDistance(unittest.TestCase):

    def setUp(self):
        self.distance = luminosity.Distance(1, 0.1)

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
            negative_distance = luminosity.Distance(-1, 0.1)

class TestBolometricLuminosity(unittest.TestCase):

    def setUp(self):
        self.lbol = luminosity.BolometricLuminosity(1E42, 1E41)

    def test_luminosity_value(self):
        expected = 1E42
        result = self.lbol.value
        self.assertEqual(expected, result)

    def test_distance_uncertainty(self):
        expected = 1E41
        result = self.lbol.uncertainty
        self.assertEqual(expected, result)

    def test_negative_distance_raises_exception(self):
        with self.assertRaises(ValueError):
            negative_lbol = luminosity.BolometricLuminosity(-1, 0.1)

class TestBolometricFlux(unittest.TestCase):

    def setUp(self):
        self.fbol = luminosity.BolometricFlux(2, 0.2)

    def test_luminosity_value(self):
        expected = 2
        result = self.fbol.value
        self.assertEqual(expected, result)

    def test_distance_uncertainty(self):
        expected = 0.2
        result = self.fbol.uncertainty
        self.assertEqual(expected, result)

    def test_negative_distance_raises_exception(self):
        with self.assertRaises(ValueError):
            negative_fbol = luminosity.BolometricFlux(-1, 0.1)

    def test_convert_flux_value_to_luminosity_value(self):
        distance = luminosity.Distance(2, 0.2)
        expected = 4.0 * math.pi * 4 * 2
        result = self.fbol.to_lbol(distance)
        self.assertEqual(expected, result.value)

    def test_convert_flux_uncertainty_to_luminosity_uncertainty(self):
        distance = luminosity.Distance(2, 0.2)
        expected = 22.479
        result = self.fbol.to_lbol(distance)
        self.assertAlmostEqual(expected, result.uncertainty, 3)

class TestConvertFluxToLuminosity(unittest.TestCase):

    def test_convert_zero_flux(self):
        expected = 0
        fbol = luminosity.BolometricFlux(0, 0)
        distance = luminosity.Distance(1, 0.1)
        result = luminosity.convert_flux_to_luminosity(fbol, distance)
        self.assertEqual(expected, result.value)

    def test_convert_flux_value_to_luminosity_value(self):
        expected = 22.479
        fbol = luminosity.BolometricFlux(2, 0.2)
        distance = luminosity.Distance(2, 0.2)
        result = luminosity.convert_flux_to_luminosity(fbol, distance)
        self.assertAlmostEqual(expected, result.uncertainty, 3)

