import unittest
import math
from .context import superbol
from superbol import mag2flux

class TestObservation(unittest.TestCase):

    def setUp(self):
        self.U_band = mag2flux.Band('U', 3660.0, 417.5E-11)
        self.B_band = mag2flux.Band('B', 4380.0, 632.0E-11)
        self.U_magnitude = 19.356
        self.B_magnitude = 18.793
        self.uncertainty = 0.02
        
    def test_observation_converts_Johnson_U_zero_magnitude_to_flux(self):
        """Observation should convert magnitude of U = 0 to the flux conversion factor of the filter"""
        magnitude = 0
        uncertainty = 0
        my_obs = mag2flux.Observation(magnitude, uncertainty, self.U_band)
        flux = my_obs.convert_to_flux()
        self.assertEqual(flux, self.U_band.flux_conversion_factor)

    def test_observation_converts_Johsnon_B_zero_magnitude_to_flux(self):
        """Observation should convert magnitude of B = 0 to the flux conversion factor of the filter"""
        magnitude = 0
        uncertainty = 0
        my_obs = mag2flux.Observation(magnitude, uncertainty, self.B_band)
        flux = my_obs.convert_to_flux()
        self.assertEqual(flux, self.B_band.flux_conversion_factor)
        
    def test_observation_converts_Johnson_U_nonzero_magnitude_to_flux(self):
        """Observation should convert magnitude of U = 19.356 to flux"""
        expected_flux = self.U_band.flux_conversion_factor * 10**(-0.4 * self.U_magnitude)
        my_obs = mag2flux.Observation(self.U_magnitude, self.uncertainty, self.U_band)
        flux = my_obs.convert_to_flux()
        self.assertEqual(flux, expected_flux)

    def test_observation_converts_Johnson_B_nonzero_magnitude_to_flux(self):
        """Observation should convert magnitude of B = 18.793 to flux"""
        expected_flux = self.B_band.flux_conversion_factor * 10**(-0.4 * self.B_magnitude)
        my_obs = mag2flux.Observation(self.B_magnitude, self.uncertainty, self.B_band)
        flux = my_obs.convert_to_flux()
        self.assertEqual(flux, expected_flux)

    def test_observation_calculates_uncertainty_in_flux_conversion(self):
        """Observation should propagate the uncertainty in the magnitude to flux conversion"""
        expected_flux = self.U_band.flux_conversion_factor * 10**(-0.4 * self.U_magnitude)
        expected_flux_uncertainty = expected_flux * 0.4 * math.log(10) * self.uncertainty
        my_obs = mag2flux.Observation(self.U_magnitude, self.uncertainty, self.U_band)
        flux_uncertainty = my_obs.calculate_flux_uncertainty()
        self.assertEqual(expected_flux_uncertainty, flux_uncertainty)

    def test_observation_zero_uncertainty(self):
        """A magnitude with no uncertainty should produce a flux with zero uncertainty"""
        uncertainty = 0.00
        expected_flux = self.B_band.flux_conversion_factor * 10**(-0.4 * self.B_magnitude)
        expected_flux_uncertainty = 0.00
        my_obs = mag2flux.Observation(self.B_magnitude, uncertainty, self.B_band)
        flux_uncertainty = my_obs.calculate_flux_uncertainty()
        self.assertEqual(expected_flux_uncertainty, flux_uncertainty)
