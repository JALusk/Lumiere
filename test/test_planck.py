import unittest
import numpy as np
from astropy import constants as const
from astropy import units as u
from planck import planck_function

class TestPlanckFunctionExtrema(unittest.TestCase):

    def setUp(self):
        self.temperature = 5000 * u.K

        # Constants that appear in BB calculations
        self.C1 = 2.0 * const.h.cgs * const.c.cgs**2
        self.C2 = const.h.cgs * const.c.cgs / const.k_B.cgs

    def test_planck_function_expected_value(self):
        wavelength = 5000 * u.Angstrom
        expected = (self.C1 / wavelength**5) / \
               (np.exp(self.C2 / (wavelength * self.temperature)) - 1) / u.sr
        expected = expected.to(u.erg / (u.s * u.cm**2 * u.AA * u.sr))
        result = planck_function(wavelength, self.temperature)

        self.assertAlmostEqual(expected.value, result.value)


    def test_rayleigh_jeans_law_long_wavelength(self):
        # The Planck function should be within 1% of the Rayleigh-Jeans
        # approximation for lambda * T > 7.7E9 angstrom * K
        wavelength = 7.8E9 * u.K * u.AA / self.temperature
        
        # Expected value is the Rayleigh-Jeans approximation
        expected = self.C1/self.C2 * self.temperature / (wavelength)**(4)
        expected = expected.to(u.erg / (u.s * u.cm**2 * u.AA))
        result = planck_function(wavelength, self.temperature)
        
        self.assertAlmostEqual(expected.value, result.value, delta = 0.01 * expected.value)

    def test_wien_approx_short_wavelength(self):
        # The Planck function should be within 1% of the Wein approximation
        # for lambda * T < 3.0E7 angstrom * K
        wavelength = 2.9E7 * u.K * u.AA / self.temperature
        
        # Expected value is the Wein approximation
        expected = self.C1 * (wavelength)**(-5) * np.exp(-self.C2 / 
                             (wavelength * self.temperature))
        expected = expected.to(u.erg / (u.s * u.cm**2 * u.AA))
        result = planck_function(wavelength, self.temperature)
        
        self.assertAlmostEqual(expected.value, result.value, delta = 0.01 * expected.value)

    def test_rayleigh_jeans_law_short_wavelength(self):
        # The Rayleigh-Jeans law should NOT work for short wavelength
        wavelength = 2.9E7 * u.K * u.AA / self.temperature

        expected = self.C1/self.C2 * self.temperature / (wavelength)**(4)
        expected = expected.to(u.erg / (u.s * u.cm**2 * u.AA))
        result = planck_function(wavelength, self.temperature)

        self.assertNotAlmostEqual(expected.value, result.value, delta = 0.01 * expected.value)

    def test_wien_approx_long_wavelength(self):
        # The Wien approximation should NOT work for long wavelength
        wavelength = 7.8E9 * u.K * u.AA / self.temperature

        expected = self.C1 * (wavelength)**(-5) * np.exp(-self.C2 / 
                             (wavelength * self.temperature))
        expected = expected.to(u.erg / (u.s * u.cm**2 * u.AA))
        result = planck_function(wavelength, self.temperature)

        self.assertNotAlmostEqual(expected.value, result.value, delta = 0.01 * expected.value)
