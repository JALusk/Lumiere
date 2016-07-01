import unittest
import numpy as np
from astropy import constants as const
from astropy import units as u
from .context import superbol
from superbol.planck import planck_function, planck_integral

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

    def test_planck_integral_UV_cool(self):
        # Taking values from Lawson, Int. J. Engng Ed. Vol. 20, No. 6, pp. 984-990, 2004
        expected = 0.000174
        result = planck_integral(3800, 2500) / (5.67051E-5 * 2500**4) * np.pi
        pdiff = abs(expected - result.value) / np.mean([expected, result.value])
        
        self.assertFalse(pdiff > 0.01)

    def test_planck_integral_UV_warm(self):
        # Taking values from Lawson, Int. J. Engng Ed. Vol. 20, No. 6, pp. 984-990, 2004
        expected = 0.052110
        result = planck_integral(3800, 5000) / (5.67051E-5 * 5000**4) * np.pi
        pdiff = abs(expected - result.value) / np.mean([expected, result.value])
        
        self.assertFalse(pdiff > 0.01)

    def test_planck_integral_UV_hot(self):
        # Taking values from Lawson, Int. J. Engng Ed. Vol. 20, No. 6, pp. 984-990, 2004
        expected = 0.443376
        result = planck_integral(3800, 10000) / (5.67051E-5 * 10000**4) * np.pi
        pdiff = abs(expected - result.value) / np.mean([expected, result.value])
        
        self.assertFalse(pdiff > 0.01)
        
    def test_planck_integral_IR_cool(self):
        # Taking values from Lawson, Int. J. Engng Ed. Vol. 20, No. 6, pp. 984-990, 2004
        expected = 0.052110
        result = planck_integral(7600, 2500) / (5.67051E-5 * 2500**4) * np.pi
        pdiff = abs(expected - result.value) / np.mean([expected, result.value])
        
        self.assertFalse(pdiff > 0.01)

    def test_planck_integral_IR_warm(self):
         # Taking values from Lawson, Int. J. Engng Ed. Vol. 20, No. 6, pp. 984-990, 2004
        expected = 0.443376
        result = planck_integral(7600, 5000) / (5.67051E-5 * 5000**4) * np.pi
        pdiff = abs(expected - result.value) / np.mean([expected, result.value])
        
        self.assertFalse(pdiff > 0.01)
        
    def test_planck_integral_IR_hot(self):
         # Taking values from Lawson, Int. J. Engng Ed. Vol. 20, No. 6, pp. 984-990, 2004
        expected = 0.839068
        result = planck_integral(7600, 10000) / (5.67051E-5 * 10000**4) * np.pi
        pdiff = abs(expected - result.value) / np.mean([expected, result.value])
        
        self.assertFalse(pdiff > 0.01)

    def test_planck_integral_convergence_R_cool(self):
        wavelength = u.Quantity(7600, unit=u.Angstrom)
        temperature = u.Quantity(2500, unit=u.K)

        C1 = 2.0 * const.h.cgs * const.c.cgs**2
        C2 = const.h.cgs * const.c.cgs / const.k_B.cgs
    
        x = C2 / (wavelength.to(u.cm) * temperature)
        iterations = min(int(2 + 20/x.value), 512)

        series = 0.0
        for i in range(1, iterations):
            term = (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6/i**4) * np.exp(-i * x)
            series += term

        B_integral = (C1 * temperature**4 / C2**4) * series / u.sr

        iterations = iterations + 1

        series = 0.0
        for i in range(1, iterations):
            term = (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6/i**4) * np.exp(-i * x)
            series += term

        B_integral_1 = (C1 * temperature**4 / C2**4) * series / u.sr

        self.assertTrue(B_integral.value / B_integral_1.value - 1 < 1E-10)

    def test_planck_integral_convergence_R_warm(self):
        wavelength = u.Quantity(7600, unit=u.Angstrom)
        temperature = u.Quantity(5000, unit=u.K)

        C1 = 2.0 * const.h.cgs * const.c.cgs**2
        C2 = const.h.cgs * const.c.cgs / const.k_B.cgs
    
        x = C2 / (wavelength.to(u.cm) * temperature)
        iterations = min(int(2 + 20/x.value), 512)

        series = 0.0
        for i in range(1, iterations):
            term = (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6/i**4) * np.exp(-i * x)
            series += term

        B_integral = (C1 * temperature**4 / C2**4) * series / u.sr

        iterations = iterations + 1

        series = 0.0
        for i in range(1, iterations):
            term = (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6/i**4) * np.exp(-i * x)
            series += term

        B_integral_1 = (C1 * temperature**4 / C2**4) * series / u.sr

        self.assertTrue(B_integral.value / B_integral_1.value - 1 < 1E-10)
        
    def test_planck_integral_convergence_R_hot(self):
        wavelength = u.Quantity(7600, unit=u.Angstrom)
        temperature = u.Quantity(10000, unit=u.K)

        C1 = 2.0 * const.h.cgs * const.c.cgs**2
        C2 = const.h.cgs * const.c.cgs / const.k_B.cgs
    
        x = C2 / (wavelength.to(u.cm) * temperature)
        iterations = min(int(2 + 20/x.value), 512)

        series = 0.0
        for i in range(1, iterations):
            term = (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6/i**4) * np.exp(-i * x)
            series += term

        B_integral = (C1 * temperature**4 / C2**4) * series / u.sr

        iterations = iterations + 1

        series = 0.0
        for i in range(1, iterations):
            term = (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6/i**4) * np.exp(-i * x)
            series += term

        B_integral_1 = (C1 * temperature**4 / C2**4) * series / u.sr

        self.assertTrue(B_integral.value / B_integral_1.value - 1 < 1E-10)

    def test_planck_integral_convergence_U_cool(self):
        wavelength = u.Quantity(21900, unit=u.Angstrom)
        temperature = u.Quantity(2500, unit=u.K)

        C1 = 2.0 * const.h.cgs * const.c.cgs**2
        C2 = const.h.cgs * const.c.cgs / const.k_B.cgs
    
        x = C2 / (wavelength.to(u.cm) * temperature)
        iterations = min(int(2 + 20/x.value), 512)

        series = 0.0
        for i in range(1, iterations):
            term = (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6/i**4) * np.exp(-i * x)
            series += term

        B_integral = (C1 * temperature**4 / C2**4) * series / u.sr

        iterations = iterations + 1

        series = 0.0
        for i in range(1, iterations):
            term = (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6/i**4) * np.exp(-i * x)
            series += term

        B_integral_1 = (C1 * temperature**4 / C2**4) * series / u.sr

        self.assertTrue(B_integral.value / B_integral_1.value - 1 < 1E-10)

    def test_planck_integral_convergence_K_warm(self):
        wavelength = u.Quantity(21900, unit=u.Angstrom)
        temperature = u.Quantity(5000, unit=u.K)

        C1 = 2.0 * const.h.cgs * const.c.cgs**2
        C2 = const.h.cgs * const.c.cgs / const.k_B.cgs
    
        x = C2 / (wavelength.to(u.cm) * temperature)
        iterations = min(int(2 + 20/x.value), 512)

        series = 0.0
        for i in range(1, iterations):
            term = (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6/i**4) * np.exp(-i * x)
            series += term

        B_integral = (C1 * temperature**4 / C2**4) * series / u.sr

        iterations = iterations + 1

        series = 0.0
        for i in range(1, iterations):
            term = (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6/i**4) * np.exp(-i * x)
            series += term

        B_integral_1 = (C1 * temperature**4 / C2**4) * series / u.sr

        self.assertTrue(B_integral.value / B_integral_1.value - 1 < 1E-10)

    def test_planck_integral_convergence_K_hot(self):
        wavelength = u.Quantity(21900, unit=u.Angstrom)
        temperature = u.Quantity(5000, unit=u.K)

        C1 = 2.0 * const.h.cgs * const.c.cgs**2
        C2 = const.h.cgs * const.c.cgs / const.k_B.cgs
    
        x = C2 / (wavelength.to(u.cm) * temperature)
        iterations = min(int(2 + 20/x.value), 512)

        series = 0.0
        for i in range(1, iterations):
            term = (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6/i**4) * np.exp(-i * x)
            series += term

        B_integral = (C1 * temperature**4 / C2**4) * series / u.sr

        iterations = iterations + 1

        series = 0.0
        for i in range(1, iterations):
            term = (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6/i**4) * np.exp(-i * x)
            series += term

        B_integral_1 = (C1 * temperature**4 / C2**4) * series / u.sr

        self.assertTrue(B_integral.value / B_integral_1.value - 1 < 1E-10)        
