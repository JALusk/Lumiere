import unittest

from unittest.mock import Mock
from unittest.mock import mock_open
from unittest.mock import patch

from .context import superbol
from superbol import read_osc
from superbol import mag2flux

file_content = """{"B":{"effective_wavelength":   4380.0,
                        "flux_conversion_factor": 632.0E-11},
                   "U":{"effective_wavelength":   3660.0,
                       "flux_conversion_factor":  417.5E-11}}"""

@patch('builtins.open', mock_open(read_data = file_content), create=True)
class TestGetObservedMagnitudes(unittest.TestCase):

    def setUp(self):
        self.magnitude = 18.78
        self.uncertainty = 0.06
        self.band_name = "B"

        self.osc_photometry_dict = {'magnitude': self.magnitude, 
                                    'e_magnitude': self.uncertainty, 
                                    'band': self.band_name}
    
    def test_no_magnitude(self):
        with self.assertRaises(read_osc.NoMagnitude):
            read_osc.get_observed_magnitude({'band': 'U'})

    def test_no_band_name_given(self):
        with self.assertRaises(read_osc.NoBandNameGiven):
            read_osc.get_observed_magnitude({'magnitude' : 0})

    def test_returns_observed_magnitude_instance(self):
        result = read_osc.get_observed_magnitude(self.osc_photometry_dict)
        self.assertIsInstance(result, mag2flux.ObservedMagnitude)

    def test_returns_correct_magnitude(self):
        result = read_osc.get_observed_magnitude(self.osc_photometry_dict)
        self.assertEqual(result.magnitude, self.magnitude)

    def test_returns_correct_uncertainty(self):
        result = read_osc.get_observed_magnitude(self.osc_photometry_dict)
        self.assertEqual(result.uncertainty, self.uncertainty)

    def test_returns_correct_band(self):
        result = read_osc.get_observed_magnitude(self.osc_photometry_dict)
        self.assertEqual(result.band.name, self.band_name)

@patch('builtins.open', mock_open(read_data = file_content), create=True)
class TestGetBand(unittest.TestCase):

    def test_returns_band_instance(self):
        result = read_osc.get_band('B')
        self.assertIsInstance(result, mag2flux.Band)

    def test_returns_B_band(self):
        band_name = 'B'
        band_effective_wavelength = 4380.0
        band_flux_conversion_factor = 632.0E-11
        result = read_osc.get_band(band_name)
        self.assertEqual(band_name, result.name)
        self.assertEqual(band_effective_wavelength, 
                         result.effective_wavelength)
        self.assertEqual(band_flux_conversion_factor, 
                         result.flux_conversion_factor)

    def test_returns_U_band(self):
        band_name = 'U'
        band_effective_wavelength = 3660.0
        band_flux_conversion_factor = 417.5E-11
        result = read_osc.get_band(band_name)
        self.assertEqual(band_name, result.name)
        self.assertEqual(band_effective_wavelength, 
                         result.effective_wavelength)
        self.assertEqual(band_flux_conversion_factor, 
                         result.flux_conversion_factor)

@patch('builtins.open', mock_open(read_data = file_content), create=True)
class TestRetrieveBandDict(unittest.TestCase):

    def test_nonexistent_band(self):
        with self.assertRaises(read_osc.NoBandFound):
            read_osc.get_band('O')

    def test_retrieves_band_data(self):
        result = read_osc.retrieve_band_dict('B')
        self.assertEqual(result["effective_wavelength"], 4380.0)
        self.assertEqual(result["flux_conversion_factor"], 632.0E-11)
