import unittest
import numpy as np
from astropy import units as u
from scipy.optimize import curve_fit
from .context import superbol
from superbol.planck import planck_function
from superbol.fit_blackbody import *

class TestFitBlackbody(unittest.TestCase):

    def setUp(self):
        self.wavelength = 5000. * u.Angstrom
        self.temperature = 13000. * u.K
        self.angular_radius = 0.2e-10
        self.flux_array = u.Quantity([5.87565760e-12, 3.79417256e-11, 
                                      6.40953095e-11, 5.90280490e-11, 
                                      3.41930932e-11], 
                                      u.erg / (u.Angstrom * u.cm**2 * u.s))
        self.flux_uncertainties = u.Quantity([1.0e-13, 1.0e-12, 
                                      1.0e-12, 1.0e-12, 
                                      1.0e-12], 
                                      u.erg / (u.Angstrom * u.cm**2 * u.s))

        self.eff_wl_array = u.Quantity([3660., 4380., 5450., 6410., 7980.],
                                        u.Angstrom)

    def test_bb_flux_returns_expected_flux(self):
        expected = (np.pi * u.sr) \
                   * planck_function(self.wavelength, self.temperature) \
                   * (self.angular_radius)**2
        result = bb_flux(self.wavelength, self.temperature, 
                         self.angular_radius)

        self.assertEqual(expected, result)
    
    def test_bb_flux_nounits_returns_expected_flux(self):
        expected = (np.pi * u.sr) \
                   * planck_function(self.wavelength, self.temperature) \
                   * (self.angular_radius)**2
        result = bb_flux_nounits(self.wavelength, self.temperature, 
                         self.angular_radius)

        self.assertEqual(expected.value, result)

    def test_bb_fit_parameters_returns_expected_parameters(self):
        popt, pcov = curve_fit(bb_flux_nounits, self.eff_wl_array.value, self.flux_array.value,
                               p0=[5000, 1.0e-10], sigma=self.flux_uncertainties.value,
                               absolute_sigma=True)
        expected_temp = popt[0]
        expected_radius = popt[1]
        expected_perr = np.sqrt(np.diag(pcov))
        result_temp, result_radius, result_perr = bb_fit_parameters(self.eff_wl_array.value,
                                                       self.flux_array.value, self.flux_uncertainties.value)
        self.assertEqual((expected_temp, expected_radius, expected_perr[0], expected_perr[1]), 
                         (result_temp, result_radius, result_perr[0], result_perr[1]))

class TestFitBlackbodyToBlackbody(unittest.TestCase):
    """Test BB fitting with actual BB flux and wavelength data"""

    def setUp(self):
        # UBVRIJHK effective wavelengths from Bessel et. al. (1998)
        self.wavelengths = [3660.0, 4380.0, 5450.0, 6410.0, 7980.0, 16300.0, 21900.0]

        # Cool, warm, and hot model temperatures
        self.cool_temp = 2500.0
        self.warm_temp = 5100.0
        self.hot_temp = 10000.0

        # Cool, warm, and hot BB fluxes at self.wavelengths
        self.cool_fluxes = [3.3786134983627783e-19, 1.825359587379945e-18, 8.074018141103508e-18, 1.7444713955287287e-17, 3.414668965720534e-17, 3.923967095683563e-17, 2.3130894084061456e-17]
        self.warm_fluxes = [1.0241915540933445e-15, 1.4831860746515914e-15, 1.7682275165286134e-15, 1.717175744676841e-15, 1.3887713520412634e-15, 2.8004081263434603e-16, 1.1313249474986263e-16] 
        self.hot_fluxes = [4.5613153169705906e-14, 3.611905858092645e-14, 2.3921603708586965e-14, 1.639405566514862e-14, 9.126671293767318e-15, 9.177243109359599e-16, 3.1983381109924263e-16]

        # Cool, warm, and hot BB flux errors (1%)
        self.cool_errors = [3.3786134983627783e-21, 1.825359587379945e-20, 8.074018141103508e-20, 1.7444713955287287e-19, 3.414668965720534e-19, 3.923967095683563e-19, 2.3130894084061456e-19]
        self.warm_errors = [4.5613153169705906e-16, 3.611905858092645e-16, 2.3921603708586965e-16, 1.639405566514862e-16, 9.126671293767318e-17, 9.177243109359599e-18, 3.1983381109924263e-18]
        self.hot_errors = [1.0241915540933445e-17, 1.4831860746515914e-17, 1.7682275165286134e-17, 1.717175744676841e-17, 1.3887713520412634e-17, 2.8004081263434603e-18, 1.1313249474986263e-18] 


        # Supernova-ish angular radiua
        self.theta = 2.0e-11
    
    def test_fit_blackbody_cool_temp(self):
        result_T, result_theta, result_perr = bb_fit_parameters(self.wavelengths, self.cool_fluxes, self.cool_errors)

        self.assertAlmostEqual(self.cool_temp, result_T)

    def test_fit_blackbody_cool_theta(self):
        result_T, result_theta, result_perr = bb_fit_parameters(self.wavelengths, self.cool_fluxes, self.cool_errors)

        self.assertAlmostEqual(self.theta, result_theta)

    def test_fit_blackbody_warm_temp(self):
        result_T, result_theta, result_perr = bb_fit_parameters(self.wavelengths, self.warm_fluxes, self.warm_errors)

        self.assertAlmostEqual(self.warm_temp, result_T)

    def test_fit_blackbody_warm_theta(self):
        result_T, result_theta, result_perr = bb_fit_parameters(self.wavelengths, self.warm_fluxes, self.warm_errors)

        self.assertAlmostEqual(self.theta, result_theta)

    def test_fit_blackbody_hot_temp(self):
        result_T, result_theta, result_perr = bb_fit_parameters(self.wavelengths, self.hot_fluxes, self.hot_errors)

        self.assertAlmostEqual(self.hot_temp, result_T)

    def test_fit_blackbody_hot_theta(self):
        result_T, result_theta, result_perr = bb_fit_parameters(self.wavelengths, self.hot_fluxes, self.hot_errors)

        self.assertAlmostEqual(self.theta, result_theta)
