import unittest
import bc_polynomial as bc_polynomial
import constants as constants

class TestSetConstants(unittest.TestCase):

    def test_set_coefficients_to_BminusV(self):
       expected = constants.coeff_BminusV
       result = bc_polynomial.set_constants("BminusV")[0]
       self.assertEqual(expected, result)
   
    def test_set_coefficients_to_VminusI(self):
       expected = constants.coeff_VminusI
       result = bc_polynomial.set_constants("VminusI")[0]
       self.assertEqual(expected, result)

    def test_set_coefficients_to_BminusI(self):
       expected = constants.coeff_BminusI
       result = bc_polynomial.set_constants("BminusI")[0]
       self.assertEqual(expected, result)

    def test_set_range_min_to_BminusV(self):
       expected = constants.min_BminusV
       result = bc_polynomial.set_constants("BminusV")[1]
       self.assertEqual(expected, result)

    def test_set_range_min_to_VminusI(self):
       expected = constants.min_VminusI
       result = bc_polynomial.set_constants("VminusI")[1]
       self.assertEqual(expected, result)

    def test_set_range_min_to_BminusI(self):
       expected = constants.min_BminusI
       result = bc_polynomial.set_constants("BminusI")[1]
       self.assertEqual(expected, result)

    def test_set_range_max_to_BminusV(self):
       expected = constants.max_BminusV
       result = bc_polynomial.set_constants("BminusV")[2]
       self.assertEqual(expected, result)

    def test_set_range_max_to_VminusI(self):
       expected = constants.max_VminusI
       result = bc_polynomial.set_constants("VminusI")[2]
       self.assertEqual(expected, result)

    def test_set_range_max_to_BminusI(self):
       expected = constants.max_BminusI
       result = bc_polynomial.set_constants("BminusI")[2]
       self.assertEqual(expected, result)

    def test_set_uncertainties_to_BminusV(self):
       expected = constants.rms_err_BminusV
       result = bc_polynomial.set_constants("BminusV")[3]
       self.assertEqual(expected, result)

    def test_set_uncertainties_to_VminusI(self):
       expected = constants.rms_err_VminusI
       result = bc_polynomial.set_constants("VminusI")[3]
       self.assertEqual(expected, result)

    def test_set_uncertainties_to_BminusI(self):
       expected = constants.rms_err_BminusI
       result = bc_polynomial.set_constants("BminusI")[3]
       self.assertEqual(expected, result)

    def test_set_constants_bad_argument_type(self):
       self.assertRaises(TypeError, bc_polynomial.set_constants, 2)

    def test_set_constants_bad_argument_value(self):
        self.assertRaises(ValueError, bc_polynomial.set_constants, 'Hello')

class TestValidityCheck(unittest.TestCase):
    
    def test_color_in_valid_range(self):
        color_value = 0.5
        range_min = 0.0
        range_max = 1.0
        self.assertTrue(bc_polynomial.valid_color(color_value,
                                                  range_min,
                                                  range_max))

    def test_color_not_in_valid_range(self):
        color_value = 2.0
        range_min = 0.0
        range_max = 1.0
        self.assertFalse(bc_polynomial.valid_color(color_value,
                                                   range_min,
                                                   range_max))
 
class TestCalculatePolynomialTerm(unittest.TestCase):

    def setUp(self):
       self.color_value = 0.5
       self.coefficient = 3.2

    def test_calculate_polynomial_term(self):
        order = 5
        expected = self.coefficient * self.color_value**(order)
        result = bc_polynomial.calculate_polynomial_term(self.coefficient, 
                                                         self.color_value,
                                                         order)
        self.assertAlmostEqual(expected, result)

    def test_calculate_polynomial_term_with_non_integer_order(self):
        order = 3.2
        self.assertRaises(TypeError, bc_polynomial.calculate_polynomial_term,
                          self.coefficient, self.color_value, order)

class TestCalculatePolynomial(unittest.TestCase):

    def setUp(self):
        self.coefficients = [1.2, 4.6, 2.5, 633.3, 34.3]
        self.variable = 3.2
    
    def test_calculate_polynomial(self):
        expected = 24390.110080000002
        result = bc_polynomial.calculate_polynomial(self.coefficients,
                                                    self.variable)
        self.assertEqual(expected, result)

class TestCalculateDerivativeTerm(unittest.TestCase):

    def setUp(self):
       self.color_value = 0.5
       self.coefficient = 3.2

    def test_calculate_polynomial_derivative_term(self):
        order = 5
        expected = order * self.coefficient * self.color_value**(order - 1)
        result = bc_polynomial.calculate_polynomial_derivative_term(
            self.coefficient, self.color_value, order)
        self.assertEqual(expected, result)

    def test_calculate_term_with_non_integer_order(self):
        order = 3.2
        self.assertRaises(TypeError, 
            bc_polynomial.calculate_polynomial_derivative_term,
            self.coefficient, self.color_value, order)

class TestQuadratureSum(unittest.TestCase):
    
    def test_quadrature_sum(self):
        expected = 5.457105459856901
        result = bc_polynomial.quadrature_sum(1.3, 5.3)
        self.assertAlmostEqual(expected,result)

class TestBolometricUncertainty(unittest.TestCase):
    
    def setUp(self):
        self.color_value = 0.5
        self.color_err = 0.04
        self.color_type = "VminusI"

    def test_bolometric_correction_err(self):
        expected = 0.11229906678151871
        result = bc_polynomial.calc_bolometric_correction_err(
                                                   self.color_value,
                                                   self.color_err,
                                                   self.color_type)
        self.assertAlmostEqual(expected, result)

class TestBolometricCorrection(unittest.TestCase):
   
    def setUp(self):
        self.color_value = 0.422
        self.color_err = 0.04
        self.color_type = "BminusV"
        
    def test_bolometric_correction(self):
        expected = -0.03984465224174367
        result = bc_polynomial.calc_bolometric_correction(self.color_value,
                                                         self.color_err,
                                                         self.color_type)[0]
        self.assertAlmostEqual(expected, result)

    def test_bolometric_correction_bad_clor_value(self):
        expected = -999
        bad_color_value = 128.54
        result = bc_polynomial.calc_bolometric_correction(bad_color_value,
                                                         self.color_err,
                                                         self.color_type)[0]
        self.assertEqual(expected, result)

if __name__ == '__main__':
    unittest.main()
