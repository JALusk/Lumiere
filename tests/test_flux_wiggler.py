import unittest

from .context import superbol
from superbol import flux_wiggler

class TestFluxWiggler(unittest.TestCase):
    def setUp(self):

    def test_wiggle_fluxes_returns_an_sed(self):
        result = flux_wiggler.wiggle_fluxes(self.sed1)
        self.assertIsInstance(result[0], mag2flux.MonochromaticFlux)

    def test_wiggle_fluxes_returns_flux_between_upper_and_lower_limit(self):
        result = flux_wiggler.wiggle_fluxes(self.sed1)
        self.assertGreaterEqual(result[0].flux, self.flux01.flux - self.flux01.flux_uncertainty)


