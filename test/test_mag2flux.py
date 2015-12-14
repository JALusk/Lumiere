import unittest
import numpy as np
from astropy import constants as const
from astropy import units as u
from mag2flux import *
from yaml import load

class TestMag2Flux(unittest.TestCase):

    def setUp(self):
        self.filter_band = "V"
        self.magnitude = 8.8
        self.uncertainty = 0.02
        self.effective_wl =  5450.0 * u.AA
        self.flux_at_zero_mag =  3.631E-9 * (u.erg / (u.s * u.cm**2 * u.AA))
    
    def test_mag2flux_converts_mag_to_correct_flux(self):
        expected = self.flux_at_zero_mag * 10**(-0.4 * self.magnitude)
        result_flux, result_uncertainty = mag2flux(self.magnitude, 
                                                              self.uncertainty,
                                                              self.effective_wl,
                                                              self.flux_at_zero_mag)

        self.assertEqual(expected.value, result_flux)

    def test_mag2flux_converts_mag_to_correct_flux_uncertainty(self):
        expected = np.sqrt((self.flux_at_zero_mag * -0.4 * np.log(10) * 10**(-0.4 * self.magnitude) * self.uncertainty)**2)
        result_flux, result_uncertainty = mag2flux(self.magnitude, 
                                                              self.uncertainty,
                                                              self.effective_wl,
                                                              self.flux_at_zero_mag)

        self.assertAlmostEqual(expected.value, result_uncertainty)


    def test_flux_at_mag_zero(self):
        mag = 0.0
        expected = self.flux_at_zero_mag
        result_flux, result_uncertainty = mag2flux(0.0,
                                                              self.uncertainty,
                                                              self.effective_wl,
                                                              self.flux_at_zero_mag)
        
        self.assertEqual(expected.value, result_flux)
