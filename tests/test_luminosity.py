import unittest
from .context import snobol
from snobol.bc_polynomial import calc_bolometric_correction as bc
import snobol.luminosity as luminosity
import snobol.constants as constants
import math


class TestLogLbol(unittest.TestCase):
    
    def setUp(self):
        self.color_value = 0.5
        self.color_err = 0.04
        self.color_type = "BminusV"
        self.v_magnitude = 16.59
        self.v_magnitude_err = 0.02
        self.distance = 1.54E23
        self.distance_err = 0.308E23
    
    def test_Fbol_calculation(self):
        expected = 10**(-0.4 * (bc(self.color_value, self.color_err,
                              self.color_type)[0]
                          + self.v_magnitude + constants.mbol_zeropoint))
        result = luminosity.calc_Fbol(self.color_value, self.color_err,
                                          self.color_type,
                                          self.v_magnitude,
                                          self.v_magnitude_err)[0]
        self.assertEqual(expected, result)

    def test_Fbol_uncertainty(self):
        Fbol = luminosity.calc_Fbol(self.color_value, self.color_err,
                                          self.color_type,
                                          self.v_magnitude,
                                          self.v_magnitude_err)[0]
        bolometric_correction, bc_err = bc(self.color_value, self.color_err, 
                                           self.color_type)
        expected = (math.sqrt(2) * 0.4 * math.log(10) * Fbol * 
                            math.sqrt(bc_err**2 + self.v_magnitude_err**2))
        result = luminosity.calc_Fbol(self.color_value, self.color_err,
                                          self.color_type,
                                          self.v_magnitude,
                                          self.v_magnitude_err)[1]
        self.assertAlmostEqual(expected, result)

    def test_Fbol_is_bad_if_bc_is_bad(self):
        color_value = 123.0
        expected = -999
        result = luminosity.calc_Fbol(color_value, self.color_err, 
                                          self.color_type,
                                          self.v_magnitude,
                                          self.v_magnitude_err)[0]
        self.assertEqual(expected, result)

    def test_Fbol_uncertainty_is_bad_if_bc_is_bad(self):
        color_value = 123.0
        expected = -999
        result = luminosity.calc_Fbol(color_value, self.color_err, 
                                          self.color_type,
                                          self.v_magnitude,
                                          self.v_magnitude_err)[1]
        self.assertEqual(expected, result)

    def test_4piDsquared_calculation(self):
        expected = 4.0 * math.pi * self.distance**2.0
        result = luminosity.calc_4piDsquared(self.distance, 
                                                 self.distance_err)[0]
        self.assertEqual(expected, result)

    def test_4piDsquared_uncertainty(self):
        expected = 8.0 * math.pi * self.distance * self.distance_err
        result = luminosity.calc_4piDsquared(self.distance,
                                                 self.distance_err)[1]
        self.assertAlmostEqual(expected, result)

    def test_Lbol(self):
        expected = luminosity.calc_Fbol(self.color_value, self.color_err, 
                                          self.color_type,
                                          self.v_magnitude,
                                          self.v_magnitude_err)[0] \
                   * luminosity.calc_4piDsquared(self.distance, 
                                               self.distance_err)[0]
        result = luminosity.calc_Lbol(self.color_value, self.color_err,
                                          self.color_type,
                                          self.v_magnitude, 
                                          self.v_magnitude_err,
                                          self.distance,
                                          self.distance_err)[0]
        self.assertAlmostEqual(expected, result)
 
    def test_Lbol_uncertainty(self):
        Fbol, Fbol_err = luminosity.calc_Fbol(self.color_value, 
                                              self.color_err, 
                                              self.color_type,
                                              self.v_magnitude,
                                              self.v_magnitude_err)
        dist, dist_err = luminosity.calc_4piDsquared(self.distance,
                                                 self.distance_err)


        expected = math.sqrt((dist * Fbol_err)**2 + (Fbol * dist_err)**2)
        result = luminosity.calc_Lbol(self.color_value, self.color_err,
                                          self.color_type,
                                          self.v_magnitude, 
                                          self.v_magnitude_err,
                                          self.distance,
                                          self.distance_err)[1]
        self.assertEqual(expected, result)
  
    def test_Lbol_is_bad_if_bc_is_bad(self):
        color_value = 123.0
        expected = -999
        result = luminosity.calc_Lbol(color_value, self.color_err,
                                          self.color_type, 
                                          self.v_magnitude, 
                                          self.v_magnitude_err,
                                          self.distance,
                                          self.distance_err)[0]
        self.assertEqual(expected, result)

    def test_Lbol_uncertainty_is_bad_if_bc_is_bad(self):
        color_value = 123.0
        expected = -999
        result = luminosity.calc_Lbol(color_value, self.color_err,
                                          self.color_type, 
                                          self.v_magnitude, 
                                          self.v_magnitude_err,
                                          self.distance,
                                          self.distance_err)[1]
        self.assertEqual(expected, result)

if __name__ == '__main__':
    unittest.main()
