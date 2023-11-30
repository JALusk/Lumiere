import unittest

from .context import superbol
from superbol import flux_wiggler
from superbol import mag2flux

class TestWiggleFluxes(unittest.TestCase):

    def setUp(self):
        #def setUp runs before each of the subsequent functions when test_flux_wiggler.py runs
        self.flux01 = mag2flux.MonochromaticFlux(100, 2, 1, 0)
        self.flux02 = mag2flux.MonochromaticFlux(200, 2, 2, 0)
        self.flux03 = mag2flux.MonochromaticFlux(150, 2, 3, 0)
        self.sed01 = [self.flux01] # sed as a list of fluxes

    def test_wiggle_fluxes_returns_an_sed(self):
        result = flux_wiggler.wiggle_fluxes(self.sed01) #Presses button in flux_wiggler to wiggle fluxes on test sed
        self.assertIsInstance(result[0], mag2flux.MonochromaticFlux) 
        # GLOSS Is result an object of the type mag2flux.MonochromaticFlux?

    def test_wiggle_fluxes_returns_a_flux_above_lower_limit(self):
        print("SED this test starts with: ")
        for flux in self.sed01:
            print(flux)
        result = flux_wiggler.wiggle_fluxes(self.sed01)
        self.assertGreaterEqual(result[0].flux, self.flux01.flux - self.flux01.flux_uncertainty) 
        # GLOSS Is the flux-first item.POSS-result object.POSS greater than or equal to the flux-test flux object.POSS-self.POSS minus the uncertainty-test flux object.POSS-self.POSS?

    def test_wiggle_fluxes_returns_a_flux_below_upper_limit(self):
        print("SED this test starts with: ")
        for flux in self.sed01:
            print(flux)
        result = flux_wiggler.wiggle_fluxes(self.sed01)
        self.assertLessEqual(result[0].flux, self.flux01.flux + self.flux01.flux_uncertainty)

    def test_wiggle_fluxes_returns_a_flux_with_the_same_uncertainty(self):
        result = flux_wiggler.wiggle_fluxes(self.sed01)
        self.assertEqual(result[0].flux_uncertainty, self.flux01.flux_uncertainty)

    def test_wiggle_fluxes_returns_a_flux_with_the_same_wavelength(self):
        result = flux_wiggler.wiggle_fluxes(self.sed01)
        self.assertEqual(result[0].wavelength, self.flux01.wavelength)

    def test_wiggle_fluxes_returns_a_flux_with_the_same_time(self):
        result = flux_wiggler.wiggle_fluxes(self.sed01)
        self.assertEqual(result[0].time, self.flux01.time)
