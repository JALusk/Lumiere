import math
import unittest

from unittest.mock import Mock
from unittest.mock import mock_open
from unittest.mock import patch

from .context import superbol
from superbol import lbc
from superbol import mag2flux
from superbol import lum


# TODO Turn this into a json file
bc_color_json_data = """{
    "H01":{
        "properties":{
            "source":{
                "bibcode":"2001PhDT.......173H",
                "reference":"Hamuy (2001)"
            },
            "ZP": -10.88802466
        },
        "B-V":{
            "range_min": -0.2,
            "range_max": 1.6,
            "coefficients": [0.199215,
                             1.654947,
                             -6.576745,
                             18.46060,
                             -25.27718,
                             15.98919,
                             -3.783559],
            "rms": 0.113
        },
        "V-I":{
            "range_min": -0.2,
            "range_max": 1.5,
            "coefficients": [-0.017371,
                             2.232705,
                             -1.246158,
                             -1.412987,
                             0.862287],
            "rms": 0.109
        }
    },
    "BH09":{
        "properties":{
            "source":{
                "bibcode":"2009ApJ...701..200B",
                "reference":"Bersten & Hamuy (2009)"
            },
            "ZP": -11.64
        },
        "B-V":{
            "range_min": -0.2,
            "range_max": 1.65,
            "coefficients": [-0.823,
                             5.027,
                             -13.409,
                             20.133,
                             -18.096,
                             9.084,
                             -1.950],
            "rms": 0.113
        },
        "V-I":{
            "range_min": -0.1,
            "range_max": 1.0,
            "coefficients": [-1.355,
                             6.262,
                             -2.676,
                             -22.973,
                             35.542,
                             -15.340],
            "rms": 0.109
        },
        "B-I":{
            "range_min": -0.4,
            "range_max": 3.0,
            "coefficients": [-1.096,
                             3.038,
                             -2.246,
                             -0.497,
                             0.7078,
                             0.576,
                             -0.713,
                             0.239,
                             -0.027],
            "rms": 0.109
        }
    }
}"""

class TestBolometricCorrection(unittest.TestCase):

    def setUp(self):
        self.BC = lbc.BolometricCorrection(3.0, 0.1)
        self.B_band = mag2flux.Band('B', 'CTIO B', 4380.0, 632.0E-11)
        self.obs = mag2flux.ObservedMagnitude(2.0, 0.1, self.B_band, 1234.5) 

    def test_apply_bolometric_correction(self):
        expected = lbc.BolometricMagnitude(5.0, 0.1, 1234.5)
        result = lbc.apply_bolometric_correction(self.BC, self.obs)
        self.assertEqual(expected.value, result.value)

class TestRetrieveBolometricCorrectionData(unittest.TestCase):

    @patch('builtins.open', mock_open(read_data = bc_color_json_data), create=True)
    def setUp(self):
        self.U_band = mag2flux.Band('U', 'CTIO U', 3660.0, 417.5E-11)
        self.B_band = mag2flux.Band('B', 'CTIO B', 4380.0, 632.0E-11)
        self.V_band = mag2flux.Band('V', 'CTIO V', 5450.0, 363.1E-11)
        self.method_name = 'BH09'
        self.method_data = lbc.retrieve_bc_method_data(self.method_name)

    @patch('builtins.open', mock_open(read_data = bc_color_json_data), create=True)
    def test_invalid_method_name(self):
        method = 'TEST'
        with self.assertRaises(lbc.InvalidBCMethod):
            result = lbc.retrieve_bc_method_data(method)

    def test_retrieve_bc_color_coefficients(self):
        expected = [-0.823, 5.027, -13.409, 20.133, -18.096, 9.084, -1.950]
        result = lbc.get_bc_method_coefficients(self.method_data, self.B_band, self.V_band)
        self.assertEqual(expected, result)

    def test_retrieve_zeropoint(self):
        expected = -11.64
        result = lbc.get_bc_method_zeropoint(self.method_data)
        self.assertEqual(expected, result)

    def test_retrieve_bc_color_coefficients_raises_invalid_filter_combination(self):
        with self.assertRaises(lbc.InvalidFilterCombination):
            result = lbc.get_bc_method_coefficients(self.method_data, self.U_band, self.B_band)

    def test_retrieve_bc_color_range(self):
        expected_min = -0.2
        expected_max = 1.65
        result = lbc.get_bc_method_range(self.method_data, self.B_band, self.V_band)
        self.assertEqual(expected_min, result[0])
        self.assertEqual(expected_max, result[1])

    def test_retrieve_bc_color_data_rms(self):
        expected = 0.113
        result = lbc.get_bc_method_rms(self.method_data, self.B_band, self.V_band)
        self.assertEqual(expected, result)

@patch('builtins.open', mock_open(read_data = bc_color_json_data), create=True)
class TestComputeBolometricCorrectionPolynomial(unittest.TestCase):

    def setUp(self):
        self.B_band = mag2flux.Band('B', 'CTIO B', 4380.0, 632.0E-11)
        self.V_band = mag2flux.Band('V', 'CTIO V', 5450.0, 363.1E-11)
        self.B_obs = mag2flux.ObservedMagnitude(18.793, 0.02, self.B_band,
                                                2451663.30)
        self.V_obs = mag2flux.ObservedMagnitude(18.078, 0.015, self.V_band,
                                                2451663.30)
        self.method = "BH09"

    def test_compute_bolometric_correction(self):
        expected = -0.0170
        result = lbc.compute_bolometric_correction(self.method, self.B_obs, self.V_obs)
        self.assertAlmostEqual(expected, result.value, 3)

    def test_compute_bolometric_correction_invalid_range(self):
        with self.assertRaises(lbc.InvalidColor):
            B_obs = mag2flux.ObservedMagnitude(18.793, 0.02, self.B_band, 2451663.30)
            V_obs = mag2flux.ObservedMagnitude(0, 0.015, self.V_band, 2451663.30)
            result = lbc.compute_bolometric_correction(self.method, B_obs, V_obs)

class TestComputePolynomial(unittest.TestCase):

    def test_zero_coefficients(self):
        coefficients = [0,0,0]
        color = 1.0
        expected = 0
        result = lbc.compute_polynomial(color, coefficients)
        self.assertEqual(expected, result)

    def test_nonzero_coefficients(self):
        coefficients = [1,2,3]
        color = 1.0
        expected = 6
        result = lbc.compute_polynomial(color, coefficients)
        self.assertEqual(expected, result)

    def test_zero_color(self):
        coefficients = [1,2,3]
        color = 0
        expected = 1
        result = lbc.compute_polynomial(color, coefficients)
        self.assertEqual(expected, result)

    def test_color_neq_one(self):
        coefficients = [1,2,3]
        color = 2.0
        expected = 1 + 2*2 + 3*2**2
        result = lbc.compute_polynomial(color, coefficients)
        self.assertEqual(expected, result)

class TestComputePolynomialDerivative(unittest.TestCase):

    def test_simple_example(self):
        coefficients = [1, 2, 3, 4]
        color = 2.0
        expected = 62
        result = lbc.compute_polynomial_derivative(color, coefficients)
        self.assertEqual(expected, result)

    def test_zero_color(self):
        coefficients = [1, 2, 3, 4]
        color = 0.0
        expected = 2
        result = lbc.compute_polynomial_derivative(color, coefficients)
        self.assertEqual(expected, result)

    def test_negative_color(self):
        coefficients = [1, 2, 3, 4]
        color = -2.0
        expected = 38
        result = lbc.compute_polynomial_derivative(color, coefficients)
        self.assertEqual(expected, result)

class TestBolometricCorrectionTechniqueH01(unittest.TestCase):

    def setUp(self):
        self.B_band = mag2flux.Band('B', 'CTIO B', 4380.0, 632.0E-11)
        self.V_band = mag2flux.Band('V', 'CTIO V', 5450.0, 363.1E-11)
        self.I_band = mag2flux.Band('I', 'CTIO I', 7980.0, 112.6E-11)
        self.B_obs = mag2flux.ObservedMagnitude(17.53, 0.015, self.B_band,
                                                2451663.30)
        self.V_obs = mag2flux.ObservedMagnitude(16.217, 0.015, self.V_band,
                                                2451663.30)
        self.I_obs = mag2flux.ObservedMagnitude(15.462, 0.015, self.I_band,
                                                2451663.30)

        self.multi_band_photometryBV = [self.B_obs, self.V_obs]
        self.multi_band_photometryVI = [self.V_obs, self.I_obs]

    def test_calculate_bc_flux_BV(self):
        expected = 7.519E-12
        result = lbc.calculate_bc_flux_h01(self.multi_band_photometryBV)
        print(expected, result.value)
        self.assertAlmostEqual(expected, result.value, delta = 0.01E-12)

    def test_calculate_bc_flux_VI(self):
        expected = 8.060E-12
        result = lbc.calculate_bc_flux_h01(self.multi_band_photometryVI)
        print(expected, result.value)
        self.assertAlmostEqual(expected, result.value, delta = 0.01E-12)
