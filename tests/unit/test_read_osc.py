import unittest

from unittest.mock import Mock
from unittest.mock import mock_open
from unittest.mock import patch

from .context import superbol
from superbol import read_osc
from superbol import mag2flux

filter_data = """{"B":{"effective_wavelength":   4380.0,
                       "flux_conversion_factor": 632.0E-11,
                       "alt_name": "CTIO B"},
                  "U":{"effective_wavelength":   3660.0,
                      "flux_conversion_factor":  417.5E-11,
                       "alt_name": "CTIO U"}}"""

@patch('builtins.open', mock_open(read_data = filter_data), create=True)
class TestGetObservedMagnitudes(unittest.TestCase):

    def setUp(self):
        self.magnitude = 18.78
        self.uncertainty = 0.06
        self.band_name = "B"
        self.time = 51663.30

        self.osc_photometry_dict = {'magnitude': "18.78", 
                                    'e_magnitude': "0.06", 
                                    'band': self.band_name,
                                    'time': self.time}
    
    def test_no_magnitude(self):
        with self.assertRaises(read_osc.NoMagnitude):
            read_osc.get_observed_magnitude({'band': 'U'})

    def test_no_band_name_given(self):
        with self.assertRaises(read_osc.NoBandNameGiven):
            read_osc.get_observed_magnitude({'magnitude' : 0})

    def test_no_time_given(self):
        with self.assertRaises(read_osc.NoTimeGiven):
            read_osc.get_observed_magnitude({'magnitude' : 0, 'band': 0})

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

    def test_returns_zero_uncertainty_if_no_uncertainty_given(self):
        result = read_osc.get_observed_magnitude({'magnitude': "18.78",
                                                  'band': self.band_name,
                                                  'time': self.time})
        self.assertEqual(result.uncertainty, 0.0)

    def test_returns_correct_time(self):
        result = read_osc.get_observed_magnitude(self.osc_photometry_dict)
        self.assertEqual(result.time, self.time)

@patch('builtins.open', mock_open(read_data = filter_data), create=True)
class TestGetBand(unittest.TestCase):

    def test_returns_band_instance(self):
        result = read_osc.get_band('B')
        self.assertIsInstance(result, mag2flux.Band)

    def test_returns_B_band(self):
        band_name = 'B'
        band_alt_name = 'CTIO B'
        band_effective_wavelength = 4380.0
        band_flux_conversion_factor = 632.0E-11
        result = read_osc.get_band(band_name)
        self.assertEqual(band_name, result.name)
        self.assertEqual(band_alt_name, result.alt_name)
        self.assertEqual(band_effective_wavelength, 
                         result.effective_wavelength)
        self.assertEqual(band_flux_conversion_factor, 
                         result.flux_conversion_factor)

    def test_returns_U_band(self):
        band_name = 'U'
        band_alt_name = 'CTIO U'
        band_effective_wavelength = 3660.0
        band_flux_conversion_factor = 417.5E-11
        result = read_osc.get_band(band_name)
        self.assertEqual(band_name, result.name)
        self.assertEqual(band_alt_name, result.alt_name)
        self.assertEqual(band_effective_wavelength, 
                         result.effective_wavelength)
        self.assertEqual(band_flux_conversion_factor, 
                         result.flux_conversion_factor)

@patch('builtins.open', mock_open(read_data = filter_data), create=True)
class TestRetrieveBandDict(unittest.TestCase):

    def test_nonexistent_band(self):
        with self.assertRaises(read_osc.NoBandFound):
            read_osc.get_band('O')

    def test_retrieves_band_data(self):
        result = read_osc.retrieve_band_dict('B')
        self.assertEqual(result["alt_name"], "CTIO B")
        self.assertEqual(result["effective_wavelength"], 4380.0)
        self.assertEqual(result["flux_conversion_factor"], 632.0E-11)

osc_json_data = """{
"test":{
        "photometry":[
                      { 
                       "time":"51663.30",
                       "band":"B",
                       "e_magnitude":"0.014",
                       "magnitude":"18.793",
                       "u_time":"MJD"
                      },
                      {
                       "time":"51663.30",
                        "band":"I",
                        "e_magnitude":"0.015",
                        "magnitude":"17.553",
                        "u_time":"MJD"
                       }
                      ]
       }
}"""

@patch('builtins.open', mock_open(read_data = osc_json_data), create=True)
class TestRetrieveOSCPhotometry(unittest.TestCase):

    def test_returns_sn_photometry(self):
        result = read_osc.retrieve_osc_photometry('test')
        self.assertEqual(result[0], {"time":"51663.30",
                                     "band":"B",
                                     "e_magnitude":"0.014",
                                     "magnitude":"18.793",
                                     "u_time":"MJD"})
