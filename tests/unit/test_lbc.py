import unittest
import math

from unittest.mock import Mock
from unittest.mock import mock_open
from unittest.mock import patch

from .context import superbol
from superbol import lbc
from superbol import mag2flux
from superbol import lqbol

from superbol.lbc import BCColorRelation

bc_color_json_data = """{
    "BH09":{
            "source":{
            "bibcode":"2009ApJ...701..200B",
            "reference":"Bersten & Hamuy (2009)"
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

class TestBCColorRelation(unittest.TestCase):

    def setUp(self):
        self.coefficients = [-0.823, 5.027, -13.409, 20.133, -18.096, 9.084, -1.950]
        self.bc_color_relation = BCColorRelation()
        self.B_band = mag2flux.Band('B', '', 0, 0)
        self.V_band = mag2flux.Band('V', '', 0, 0)
        self.I_band = mag2flux.Band('I', '', 0, 0)

    def test_get_coefficients_BV(self):
        result = self.bc_color_relation.get_coefficients(self.B_band, self.V_band)
        self.assertEqual(self.coefficients, result)

    def test_get_coefficients_VI(self):
        expected = [-1.355, 6.262, -2.676, -22.973, 35.542, -15.340]
        result = self.bc_color_relation.get_coefficients(self.V_band, self.I_band)
        self.assertEqual(expected, result)

    def get_bolometric_correction(self):
        bc = self.bc_color_relation.get_bolometric_correction(self.B_obs, self.V_obs)
        expected = -0.0170
        self.assertEqual(expected, bc.value)

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

class TestBCLuminosity(unittest.TestCase):

    def setUp(self):
        # IAU 2015 Resolution B2 value (in erg/s)
        self.L0 = 3.0128E35

    def test_convert_zero_magnitude_to_bolometric_luminosity(self):
        expected = self.L0
        Mbol_zero = lbc.BolometricMagnitude(0, 0, 0)
        result = lbc.convert_Mbol_to_Lbol(Mbol_zero)
        self.assertEqual(expected, result.value)

    def test_convert_solar_magnitude_to_bolometric_luminosity(self):
        Msun = lbc.BolometricMagnitude(4.74, 0, 0)
        expected = self.L0 * 10**(-0.4 * Msun.value)
        result = lbc.convert_Mbol_to_Lbol(Msun)
        self.assertEqual(expected, result.value)

    def test_convert_magnitude_to_bolometric_luminosity_uncertainty(self):
        Mtest = lbc.BolometricMagnitude(5.0, 0.1, 0)
        expected = -0.4 * self.L0 * 10**(-0.4 * Mtest.value)

class TestDistanceModulus(unittest.TestCase):

    def setUp(self):
        self.msun = lbc.BolometricMagnitude(-26.832,0,0)
        self.dsun = lqbol.Distance(1.496E13, 0.0)

    def test_convert_solar_apparent_to_absolute_magnitude(self):
        expected = 4.74
        result = lbc.convert_apparent_to_absolute_magnitude(self.msun, self.dsun)
        self.assertAlmostEqual(expected, result.value, 3)
