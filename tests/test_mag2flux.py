import unittest
from .context import superbol
from superbol import mag2flux

class TestObservation(unittest.TestCase):

    def test_observation(self):
        my_obs = mag2flux.Observation(0)

    def test_observation_converts_Johnson_U_zero_magnitude_to_flux(self):
        """Observation should convert magnitude of U = 0 to the flux conversion factor of the filter"""
        magnitude = 0
        flux_conversion_factor = 417.5E-11
        my_obs = mag2flux.Observation(magnitude)
        flux = my_obs.convert_to_flux()
        self.assertEqual(flux, flux_conversion_factor)
        
