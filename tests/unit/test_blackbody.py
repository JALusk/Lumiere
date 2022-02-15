import unittest
import numpy as np
from astropy import units as u
from .context import superbol
from superbol import mag2flux
from superbol.planck import planck_function
from superbol.blackbody import *

# TODO May combine these into one
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

        self.blackbody_fit = BlackbodyFit()
        self.blackbody_fit.fit_to_SED(self.SED)

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
