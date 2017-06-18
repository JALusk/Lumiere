import unittest
import math
from .context import superbol
from superbol import mag2flux

class TestObservation(unittest.TestCase):

    def test_observation_converts_Johnson_U_zero_magnitude_to_flux(self):
        """Observation should convert magnitude of U = 0 to the flux conversion factor of the filter"""
        magnitude = 0
        uncertainty = 0
        flux_conversion_factor = 417.5E-11
        my_obs = mag2flux.Observation(magnitude, uncertainty, flux_conversion_factor)
        flux = my_obs.convert_to_flux()
        self.assertEqual(flux, flux_conversion_factor)

    def test_observation_converts_Johsnon_B_zero_magnitude_to_flux(self):
        """Observation should convert magnitude of B = 0 to the flux conversion factor of the filter"""
        magnitude = 0
        uncertainty = 0
        flux_conversion_factor = 632.0E-11
        my_obs = mag2flux.Observation(magnitude, uncertainty, flux_conversion_factor)
        flux = my_obs.convert_to_flux()
        self.assertEqual(flux, flux_conversion_factor)
        
    def test_observation_converts_Johnson_U_nonzero_magnitude_to_flux(self):
        """Observation should convert magnitude of U = 19.356 to flux"""
        magnitude = 19.356
        uncertainty = 0.02
        flux_conversion_factor = 417.5E-11
        expected_flux = flux_conversion_factor * 10**(-0.4 * magnitude)
        my_obs = mag2flux.Observation(magnitude, uncertainty, flux_conversion_factor)
        flux = my_obs.convert_to_flux()
        self.assertEqual(flux, expected_flux)

    def test_observation_converts_Johnson_B_nonzero_magnitude_to_flux(self):
        """Observation should convert magnitude of B = 18.793 to flux"""
        magnitude = 18.793
        uncertainty = 0.02
        flux_conversion_factor = 632.0E-11
        expected_flux = flux_conversion_factor * 10**(-0.4 * magnitude)
        my_obs = mag2flux.Observation(magnitude, uncertainty, flux_conversion_factor)
        flux = my_obs.convert_to_flux()
        self.assertEqual(flux, expected_flux)

    def test_observation_calculates_uncertainty_in_flux_conversion(self):
        """Observation should propagate the uncertainty in the magnitude to flux conversion"""
        magnitude = 19.356
        uncertainty = 0.02
        flux_conversion_factor = 417.5E-11
        expected_flux = flux_conversion_factor * 10**(-0.4 * magnitude)
        expected_flux_uncertainty = expected_flux * 0.4 * math.log(10) * uncertainty
        my_obs = mag2flux.Observation(magnitude, uncertainty, flux_conversion_factor)
        flux_uncertainty = my_obs.calculate_flux_uncertainty()
        self.assertEqual(expected_flux_uncertainty, flux_uncertainty)

    def test_observation_zero_uncertainty(self):
        """A magnitude with no uncertainty should produce a flux with zero uncertainty"""
        magnitude = 19.356
        uncertainty = 0.00
        flux_conversion_factor = 417.5E-11
        expected_flux = flux_conversion_factor * 10**(-0.4 * magnitude)
        expected_flux_uncertainty = 0.00
        my_obs = mag2flux.Observation(magnitude, uncertainty, flux_conversion_factor)
        flux_uncertainty = my_obs.calculate_flux_uncertainty()
        self.assertEqual(expected_flux_uncertainty, flux_uncertainty)
        
