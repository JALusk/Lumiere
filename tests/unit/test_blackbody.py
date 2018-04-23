import unittest
import numpy as np
from astropy import units as u
from scipy.optimize import curve_fit
from .context import superbol
from superbol import mag2flux
from superbol import sed
from superbol.planck import planck_function
from superbol.blackbody import *
from scipy import integrate

class TestFitBlackbody(unittest.TestCase):

    def setUp(self):
        self.wavelength = 5000. * u.Angstrom
        self.temperature = 13000. * u.K
        self.angular_radius = 0.2e-10
        self.flux1 = mag2flux.MonochromaticFlux(flux = 5.87565760e-12,
                                                flux_uncertainty = 1.0e-13,
                                                wavelength = 3600.0,
                                                time = 0)
        self.flux2 = mag2flux.MonochromaticFlux(flux = 3.79417256e-11,
                                                flux_uncertainty = 1.0e-12,
                                                wavelength = 4380.0,
                                                time = 0)
        self.flux3 = mag2flux.MonochromaticFlux(flux = 6.40953095e-11,
                                                flux_uncertainty = 1.0e-12,
                                                wavelength = 5450.0,
                                                time = 0)
        self.flux4 = mag2flux.MonochromaticFlux(flux = 5.90280490e-11,
                                                flux_uncertainty = 1.0e-12,
                                                wavelength = 6410.0,
                                                time = 0)
        self.flux5 = mag2flux.MonochromaticFlux(flux = 3.41930932e-11,
                                                flux_uncertainty = 1.0e-12,
                                                wavelength = 7980.0,
                                                time = 0)
        self.SED = [self.flux1, self.flux2, self.flux3, self.flux4, self.flux5]

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

#    def test_bb_fit_parameters_returns_expected_parameters(self):
#        wavelengths = [f.wavelength for f in self.SED]
#        fluxes = sed.get_flux_values(self.SED)
#        flux_uncertainties = sed.get_flux_uncertainties(self.SED)
#        popt, pcov = curve_fit(bb_flux_nounits, wavelengths, fluxes,
#                               p0=[5000, 1.0e-10], sigma=flux_uncertainties,
#                               absolute_sigma=True)
#        expected_temp = popt[0]
#        expected_radius = popt[1]
#        expected_perr = np.sqrt(np.diag(pcov))
#        result_temp, result_radius, result_perr = bb_fit_parameters(self.SED)
#        self.assertEqual((expected_temp, expected_radius, 
#                          expected_perr[0], expected_perr[1]), 
#                         (result_temp, result_radius, 
#                          result_perr[0], result_perr[1]))
#
#class TestFitBlackbodyToBlackbody(unittest.TestCase):
#    """Test BB fitting with actual BB flux and wavelength data"""
#
#    def setUp(self):
#        # UBVRIJHK effective wavelengths from Bessel et. al. (1998)
#        self.wavelengths = [3660.0, 4380.0, 5450.0, 6410.0, 7980.0, 16300.0, 21900.0]
#
#        # Cool, warm, and hot model temperatures
#        self.cool_temp = 2500.0
#        self.warm_temp = 5100.0
#        self.hot_temp = 10000.0
#
#        # Cool, warm, and hot BB fluxes at self.wavelengths
#        self.cool_fluxes = [3.3786134983627783e-19, 1.825359587379945e-18, 8.074018141103508e-18, 1.7444713955287287e-17, 3.414668965720534e-17, 3.923967095683563e-17, 2.3130894084061456e-17]
#        self.warm_fluxes = [1.0241915540933445e-15, 1.4831860746515914e-15, 1.7682275165286134e-15, 1.717175744676841e-15, 1.3887713520412634e-15, 2.8004081263434603e-16, 1.1313249474986263e-16] 
#        self.hot_fluxes = [4.5613153169705906e-14, 3.611905858092645e-14, 2.3921603708586965e-14, 1.639405566514862e-14, 9.126671293767318e-15, 9.177243109359599e-16, 3.1983381109924263e-16]
#
#        # Cool, warm, and hot BB flux errors (1%)
#        self.cool_errors = [3.3786134983627783e-21, 1.825359587379945e-20, 8.074018141103508e-20, 1.7444713955287287e-19, 3.414668965720534e-19, 3.923967095683563e-19, 2.3130894084061456e-19]
#        self.warm_errors = [4.5613153169705906e-16, 3.611905858092645e-16, 2.3921603708586965e-16, 1.639405566514862e-16, 9.126671293767318e-17, 9.177243109359599e-18, 3.1983381109924263e-18]
#        self.hot_errors = [1.0241915540933445e-17, 1.4831860746515914e-17, 1.7682275165286134e-17, 1.717175744676841e-17, 1.3887713520412634e-17, 2.8004081263434603e-18, 1.1313249474986263e-18]
#
#        # Build SEDs
#        self.cool_SED = []
#        self.warm_SED = []
#        self.hot_SED = []
#        for i in range(len(self.wavelengths)):
#            cool_flux = mag2flux.MonochromaticFlux(self.cool_fluxes[i], self.cool_errors[i], self.wavelengths[i], 0)
#            warm_flux = mag2flux.MonochromaticFlux(self.warm_fluxes[i], self.warm_errors[i], self.wavelengths[i], 0)
#            hot_flux = mag2flux.MonochromaticFlux(self.hot_fluxes[i], self.hot_errors[i], self.wavelengths[i], 0)
#            self.cool_SED.append(cool_flux)
#            self.warm_SED.append(warm_flux)
#            self.hot_SED.append(hot_flux)
#
#        # Supernova-ish angular radiua
#        self.theta = 2.0e-11
#    
#    def test_fit_blackbody_cool_temp(self):
#        result_T, result_theta, result_perr = bb_fit_parameters(self.cool_SED)
#
#        self.assertAlmostEqual(self.cool_temp, result_T, 2)
#
#    def test_fit_blackbody_cool_theta(self):
#        result_T, result_theta, result_perr = bb_fit_parameters(self.cool_SED)
#
#        self.assertAlmostEqual(self.theta, result_theta)
#
#    def test_fit_blackbody_warm_temp(self):
#        result_T, result_theta, result_perr = bb_fit_parameters(self.warm_SED)
#
#        self.assertAlmostEqual(self.warm_temp, result_T, 2)
#
#    def test_fit_blackbody_warm_theta(self):
#        result_T, result_theta, result_perr = bb_fit_parameters(self.warm_SED)
#
#        self.assertAlmostEqual(self.theta, result_theta)
#
#    def test_fit_blackbody_hot_temp(self):
#        result_T, result_theta, result_perr = bb_fit_parameters(self.hot_SED)
#
#        self.assertAlmostEqual(self.hot_temp, result_T, 2)
#
#    def test_fit_blackbody_hot_theta(self):
#        result_T, result_theta, result_perr = bb_fit_parameters(self.hot_SED)
#
#        self.assertAlmostEqual(self.theta, result_theta)
#
#class TestFitBlackbodyTemperatureToWD(unittest.TestCase):
#    """Test BB fitting with WD flux and wavelength data"""
#
#    def setUp(self):
#        # Effective wavelengths
#        self.wavelengths = [3639.3, 3660.0, 4380.0,
#                            4765.1, 5450.0, 6223.3,
#                            6410.0, 7609.2, 7980.0, 
#                            12200.0, 16300.0, 21900.0]
#
#        # Cool, warm, and hot model temperatures
#        self.cool_temp = 3500.0
#        self.warm_temp = 5000.0
#        self.hot_temp = 10000.0
#
#        # Cool, warm, and hot BB fluxes at self.wavelengths
#        self.cool_fluxes = [3.17400826e-17, 3.34700591e-17, 1.16712421e-16,
#                            1.86464465e-16, 3.02291718e-16, 4.28431996e-16,
#                            4.58203528e-16, 5.21624341e-16, 5.59678199e-16,
#                            3.91108039e-16, 2.18847839e-16, 1.00139367e-16]
#        self.warm_fluxes = [1.53889168e-15, 1.56695484e-15, 2.51140687e-15,
#                            2.94982383e-15, 3.20941774e-15, 3.32875550e-15,
#                            3.12649997e-15, 2.75776877e-15, 2.59381138e-15,
#                            1.09420008e-15, 4.95840803e-16, 1.97244473e-16]
#        self.hot_fluxes = [9.62073059e-14, 9.81423561e-14, 6.84096043e-14,
#                           6.04499687e-14, 4.30552632e-14, 3.30431846e-14,
#                           2.81746435e-14, 1.88867769e-14, 1.56869456e-14,
#                           4.02806753e-15, 1.47279487e-15, 5.04666925e-16]
#
#        # Cool, warm, and hot BB flux errors (1%)
#        self.cool_errors = [5.84673929e-19, 6.16541274e-19, 2.14992225e-18,
#                            3.43480238e-18, 5.56841923e-18, 7.89200901e-18,
#                            8.44042091e-18, 9.60867546e-18, 1.03096534e-17,
#                            7.20447633e-18, 4.03132617e-18, 1.84463531e-18]
#        self.warm_errors = [2.83474324e-17, 2.88643748e-17, 4.62618242e-17,
#                            5.43377631e-17, 5.91196596e-17, 6.13179423e-17,
#                            5.75922578e-17, 5.07999781e-17, 4.77797713e-17,
#                            2.01559104e-17, 9.13372513e-18, 3.63337746e-18]
#        self.hot_errors = [1.77220407e-15, 1.80784901e-15, 1.26015148e-15,
#                           1.11352957e-15, 7.93107259e-16, 6.08677954e-16,
#                           5.18996112e-16, 3.47907287e-16, 2.88964217e-16,
#                           7.41997460e-17, 2.71298841e-17, 9.29630830e-18]
#
#        # Build SEDs
#        self.cool_SED = []
#        self.warm_SED = []
#        self.hot_SED = []
#        for i in range(len(self.wavelengths)):
#            cool_flux = mag2flux.MonochromaticFlux(self.cool_fluxes[i], self.cool_errors[i], self.wavelengths[i], 0)
#            warm_flux = mag2flux.MonochromaticFlux(self.warm_fluxes[i], self.warm_errors[i], self.wavelengths[i], 0)
#            hot_flux = mag2flux.MonochromaticFlux(self.hot_fluxes[i], self.hot_errors[i], self.wavelengths[i], 0)
#            self.cool_SED.append(cool_flux)
#            self.warm_SED.append(warm_flux)
#            self.hot_SED.append(hot_flux)
#
#
#    def test_fit_blackbody_cool_WD(self):
#        result_T, result_theta, result_perr = bb_fit_parameters(self.cool_SED)
#
#        expected = 3300.274543
#        self.assertAlmostEqual(expected, result_T, 2)
#
#    def test_fit_blackbody_warm_WD(self):
#        result_T, result_theta, result_perr = bb_fit_parameters(self.warm_SED)
#
#        expected = 4983.36073925
#        self.assertAlmostEqual(expected, result_T, 2)
#
#    def test_fit_blackbody_hot_WD(self):
#        result_T, result_theta, result_perr = bb_fit_parameters(self.hot_SED)
#
#        expected = 11196.3769245
#        self.assertAlmostEqual(expected, result_T, 2)
#
#
#class TestFitBlackbodyTemperatureToSN(unittest.TestCase):
#
#    def setUp(self):
#        self.wavelengths = [3660., 4380.,  5450.,  6410.,  7980.]
#        self.fluxes = [1.42159582e-16, 2.09548255e-16, 3.70796730e-16, 4.53049602e-16, 3.36788260e-16]
#        self.flux_errs = [2.21846926e-17, 1.18302183e-17, 1.17834356e-17, 1.53145455e-17, 1.24389601e-17]
#
#        self.SED = []
#        for i in range(len(self.wavelengths)):
#            flux = mag2flux.MonochromaticFlux(self.fluxes[i], self.flux_errs[i], self.wavelengths[i], 0)
#            self.SED.append(flux)
#
#    def test_fit_blackbody_SN98A(self):
#        result_T, result_theta, result_perr = bb_fit_parameters(self.SED)
#
#        expected = 4331.7954035003195
#        self.assertAlmostEqual(expected, result_T, 2)
#
#class TestEmptyBlackbodyFit(unittest.TestCase):
#
#    def setUp(self):
#        self.empty_blackbody_fit = BlackbodyFit()
#
#    def test_empty_blackbody_fit_temperature(self):
#        expected = None
#        result = self.empty_blackbody_fit.temperature
#        self.assertEqual(expected, result)
#   
#    def test_empty_blackbody_fit_angular_radius(self):
#        expected = None
#        result = self.empty_blackbody_fit.angular_radius
#        self.assertEqual(expected, result)
#
#    def test_empty_blackbody_fit_temperature_err(self):
#        expected = None
#        result = self.empty_blackbody_fit.temperature_err
#        self.assertEqual(expected, result)
#
#    def test_empty_blackbody_fit_angular_radius_err(self):
#        expected = None
#        result = self.empty_blackbody_fit.angular_radius_err
#        self.assertEqual(expected, result)
#
#    def test_empty_blackbody_fit_SED(self):
#        expected = None
#        result = self.empty_blackbody_fit.SED
#        self.assertEqual(expected, result)
#
#    def test_empty_blackbody_fit_call(self):
#        expected = None
#        result = self.empty_blackbody_fit(5000)
#        self.assertEqual(expected, result)

class TestBlackbodyFit(unittest.TestCase):

    def setUp(self):
        self.wavelengths = [3660.0, 4380.0, 5450.0, 6410.0, 
                            7980.0, 16300.0, 21900.0]
        self.fluxes = [1.0241915540933445e-15, 1.4831860746515914e-15, 
                       1.7682275165286134e-15, 1.717175744676841e-15, 
                       1.3887713520412634e-15, 2.8004081263434603e-16,
                       1.1313249474986263e-16] 
        self.errors = [4.5613153169705906e-16, 3.611905858092645e-16, 
                       2.3921603708586965e-16, 1.639405566514862e-16,
                       9.126671293767318e-17, 9.177243109359599e-18, 
                       3.1983381109924263e-18]
        
        self.SED = []
        for i in range(len(self.wavelengths)):
            flux = mag2flux.MonochromaticFlux(self.fluxes[i], self.errors[i], 
                                              self.wavelengths[i], 0)
            self.SED.append(flux)

        self.blackbody_fit = BlackbodyFit()
        self.blackbody_fit.fit_to_SED(self.SED)
 
    def test_blackbody_fit_temperature_to_SED(self):
        expected = 5100.0
        result = self.blackbody_fit.temperature
        self.assertAlmostEqual(expected, result, 2)

    def test_blackbody_fit_angular_radius_to_SED(self):
        expected = 2.0E-11
        result = self.blackbody_fit.angular_radius
        self.assertAlmostEqual(expected, result, 13)

    def test_blackbody_fit_I_flux(self):
        wavelength = 7980.0
        temperature = 5100.0 * u.K
        angular_radius = 2.0E-11
        expected = 1.3887713520412634e-15
        result = self.blackbody_fit(wavelength)
        self.assertAlmostEqual(expected, result, 17)

    def test_blackbody_fit_V_flux(self):
        wavelength = 5450.0
        temperature = 5100.0 * u.K
        angular_radius = 2.0E-11
        expected = 1.7682275165286134e-15
        result = self.blackbody_fit(wavelength)
        self.assertAlmostEqual(expected, result, 17)
