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
