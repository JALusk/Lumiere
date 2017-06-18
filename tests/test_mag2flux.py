import unittest
from .context import superbol
from superbol import mag2flux

class TestObservation(unittest.TestCase):

    def test_observation_converts_Johnson_U_zero_magnitude_to_flux(self):
        """Observation should convert magnitude of U = 0 to the flux conversion factor of the filter"""
        magnitude = 0
        flux_conversion_factor = 417.5E-11
        my_obs = mag2flux.Observation(magnitude, flux_conversion_factor)
        flux = my_obs.convert_to_flux()
        self.assertEqual(flux, flux_conversion_factor)

    def test_observation_converts_Johsnon_B_zero_magnitude_to_flux(self):
        """Observation should convert magnitude of B = 0 to the flux conversion factor of the filter"""
        magnitude = 0
        flux_conversion_factor = 632.0E-11
        my_obs = mag2flux.Observation(magnitude, flux_conversion_factor)
        flux = my_obs.convert_to_flux()
        self.assertEqual(flux, flux_conversion_factor)
        
    def test_observation_converts_Johnson_U_nonzero_magnitude_to_flux(self):
        """Observation should convert magnitude of U = 19.356 to flux"""
        magnitude = 19.356
        flux_conversion_factor = 417.5E-11
        expected_flux = flux_conversion_factor * 10**(-0.4 * magnitude)
        my_obs = mag2flux.Observation(magnitude, flux_conversion_factor)
        flux = my_obs.convert_to_flux()
        self.assertEqual(flux, expected_flux)

    def test_observation_converts_Johnson_B_nonzero_magnitude_to_flux(self):
        """Observation should convert magnitude of B = 18.793 to flux"""
        magnitude = 18.793
        flux_conversion_factor = 632.0E-11
        expected_flux = flux_conversion_factor * 10**(-0.4 * magnitude)
        my_obs = mag2flux.Observation(magnitude, flux_conversion_factor)
        flux = my_obs.convert_to_flux()
        self.assertEqual(flux, expected_flux)
