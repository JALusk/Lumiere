import unittest

from .context import superbol
from superbol import laug
from superbol import mag2flux
from superbol import blackbody

class TestTrimSED(unittest.TestCase):

    def setUp(self):
        self.time = 1234.5
        self.flux1 = mag2flux.MonochromaticFlux(flux = 1,
                                                flux_uncertainty = 1,
                                                wavelength = 3660.0,
                                                time = self.time)
        self.flux2 = mag2flux.MonochromaticFlux(flux = 2,
                                                flux_uncertainty = 1,
                                                wavelength = 4380.0,
                                                time = self.time)
        self.flux3 = mag2flux.MonochromaticFlux(flux = 3,
                                                flux_uncertainty = 1,
                                                wavelength = 5450.0,
                                                time = self.time)
        self.flux4 = mag2flux.MonochromaticFlux(flux = 4,
                                                flux_uncertainty = 1,
                                                wavelength = 6410.0,
                                                time = self.time)
        self.flux5 = mag2flux.MonochromaticFlux(flux = 5,
                                                flux_uncertainty = 1,
                                                wavelength = 7980.0,
                                                time = self.time)
        self.SED = [self.flux1, self.flux2, self.flux3, self.flux4, self.flux5]

    def test_trim_SED_5000_angstroms(self):
        expected = [self.flux3, self.flux4, self.flux5]
        result = laug.trim_SED(self.SED, 5000.0)
        
        self.assertEqual(expected, result)

    def test_trim_SED_7000_angstroms(self):
        expected = [self.flux5]
        result = laug.trim_SED(self.SED, 7000.0)

        self.assertEqual(expected, result)

class TestIRCorrection(unittest.TestCase):
    pass
