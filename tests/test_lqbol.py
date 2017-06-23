import unittest
import math
from .context import superbol
from superbol import mag2flux
from superbol import lqbol

class TestSpectralEnergyDistribution(unittest.TestCase):

    def setUp(self):
        self.sed = lqbol.SpectralEnergyDistribution()

    def test_sort_dummy_flux_objects_by_wavelength(self):
        flux2 = mag2flux.MonochromaticFlux(300, 10, 2)
        flux1 = mag2flux.MonochromaticFlux(200, 30, 1)
        flux3 = mag2flux.MonochromaticFlux(100, 20, 3)
        fluxes = [flux2, flux1, flux3]
        expected = [flux1, flux2, flux3]
        result = self.sed.sort_by_wavelength(fluxes)
        self.assertEqual(expected, result)

    def test_sort_real_flux_objects_by_wavelength(self):
        flux1 = mag2flux.MonochromaticFlux(1.921E-16, 2.654E-18, 4380.0)
        flux2 = mag2flux.MonochromaticFlux(7.555E-17, 2.505E-18, 3660.0)
        flux3 = mag2flux.MonochromaticFlux(2.132E-16, 2.946E-18, 5450.0)
        fluxes = [flux1, flux2, flux3]
        expected = [flux2, flux1, flux3]
        result = self.sed.sort_by_wavelength(fluxes)
        self.assertEqual(expected, result)
