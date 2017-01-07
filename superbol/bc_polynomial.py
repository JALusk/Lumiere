import superbol.constants as constants
import math

def set_constants(color_type):
    """Sets the coefficients, validty range, and rms error of fit.

    Args:
        color_type (str): A string specifying the color combination. Must be
            "BminusV" for B-V, "VminusI" for V-I, or "BminusI" for B-I.

    Returns:
        tuple: A tuple containing the list of coefficients for the polynomial
        fit which correspond to the supplied color, the minimum value of
        the color for which the polynomial is valid, the maximum value
        of the color for which the polynomial is valid, and the rms
        error of the polynomial fit for the supplied color.

        ([coefficients], min, max, rms_error)

    Raises:
        TypeError: The argument given is not a string
        ValueError: The argument given is not one of the three valid
            strings.
    """
    if color_type == "BminusV":
        return constants.coeff_BminusV, constants.min_BminusV, \
               constants.max_BminusV, constants.rms_err_BminusV
    elif color_type == "VminusI":
        return constants.coeff_VminusI, constants.min_VminusI, \
               constants.max_VminusI, constants.rms_err_VminusI
    elif color_type == "BminusI":
        return constants.coeff_BminusI, constants.min_BminusI, \
               constants.max_BminusI, constants.rms_err_BminusI
    elif type(color_type) != str:
        raise TypeError("The argument given is not a string")
    else:
        raise ValueError("The argument given is not a valid color")

def valid_color(color_value, range_min, range_max):
    """Checks that the color value is within the range of validity.

    Each polynomial fit has a different range of validity. We need to
    make sure that the color we are feeding in is inside that range.

    Args:
        color_value (float): B-V, V-I, or B-I color of the supernova in
            magnitudes (corrected for reddening and extinction from the
            host and MWG.)
        range_min (float): Minumum value of the color range over which the fit
            is valid
        range_max (float): Maxumum value of the color range over which the fit
            is valid

    Returns:
        bool: True if the color value is within the valid range.
        False if the volor value is outside the valid range.
   """
    if range_min <= color_value <= range_max:
        return True
    else:
        return False

def calculate_polynomial_term(coefficient, variable, order):
    """Calculates a term in a polynomial.

    Args:
        coefficient (float): The coefficient to use in
            calculating the polynomial term.
        variable (float): Value to plug in for the variable in the polynomial
            term.
        order (int): Integer to use as the order of the polynomial term.

    Returns:
        float: The result of coefficient * variable**(order)

    Raises:
        TypeError: A non-integer order is given.
    """
    if type(order) != int:
        raise TypeError('Non-integer order in polynomial term')
    else:
        return coefficient * variable**(order)

def calculate_polynomial(coefficients, variable):
    """Calculates a polynomial.

    Args:
        coefficients (list): list of polynomial coefficients. The length
            of the list will be used as the order of the polynomial.
        variable (float): float to plug in for the variable in the polynomial.

    Returns:
        float: The result of summing the polynomial terms
        calculated from the coefficients and variable given.
    """
    polynomial = 0.0

    for order in range(len(coefficients)):
        polynomial += calculate_polynomial_term(coefficients[order],
                                                variable,
                                                order)

    return polynomial

def calculate_polynomial_derivative_term(coefficient, variable, order):
    """Calculates the derivative of the nth order term of a polynomial.

    Args:
        coefficient (float): The coefficient of the nth order term in the
            polynomial
        variable (float): float to plug in for the variable in the polynomial
        order (int): order of the nth order term in the polynomial (so, n.)

    Returns:
        float: The result of taking the derivative of the nth
        order term a polynomial,

        :math:`n \\cdot \\text{coefficient} \\cdot \\text{variable}^{n-1}`

        So, the edge case of taking the derivative of the zeroth-order
        term is taken care of, since you explicity multiply by the
        order of the polynomial (which is zero in the n = 0 case.)

    Raises:
        TypeError: A non-integer was passed as the order.
    """
    if type(order) != int:
        raise TypeError('Non-integer order in polynomial term')
    else:
        return order * coefficient * variable**(order - 1)

def calculate_polynomial_derivative(coefficients, variable):
    """Calculates the derivative of a polynomial.

    Args:
        coefficients (list): List of polynomial coefficients. The length
            of the list will be used as the order of the polynomial.
        variable (float): Value to plug in for the variable in the polynomial.

    Returns:
        float: The result of summing the derivatives of
        the polynomial terms calculated from the coefficients and
        variable given.
    """
    polynomial_derivative = 0.0

    for order in range(len(coefficients)):
        polynomial_derivative += calculate_polynomial_derivative_term(
            coefficients[order], variable, order)

    return polynomial_derivative

def quadrature_sum(x, y):
    """Calculate the quadrature sum of two variables x and y.

    Args:
        x (float): Variable to include in the quadrature sum
        y (float): Variable to include in the quadrature sum

    Returns:
        float: Result of calculating :math:`\\sqrt{x^2 + y^2}`
    """
    return math.sqrt(x**2 + y**2)

def calc_bolometric_correction_err(color_value, color_err, color_type):
    """Calculates the uncertainty in the bolometric correction.

    Two uncertainties are added in quadrature to get the total
    uncertainty in the bolometric correction. The first is uncertainty
    in the BC due to uncertainties in the measured color value (simple
    error propagation using a derivative.) The second is the RMS error
    inherent in the polynomial fit to the template data as reported in
    Bersten & Hamuy (2009.)

    Args:
        color_value (float): B-V, V-I, or B-I color of the supernova in
            magnitudes (corrected for reddening and extinction from
            the host and MWG.)
        color_err (float): Uncertainty in the photometric color.
        color_type (str): String signifying which color color_value
            represents.

    Returns:
        float: Uncertainty in the value of the bolometric correction
   """
    coefficients = set_constants(color_type)[0]
    rms_err = set_constants(color_type)[3]

    bc_derivative = calculate_polynomial_derivative(coefficients,
                                                     color_value)
    bc_polynomial_err = abs(bc_derivative) * color_err
    bolometric_correction_uncertainty = quadrature_sum(bc_polynomial_err,
                                                       rms_err)

    return bolometric_correction_uncertainty

def calc_bolometric_correction(color_value, color_err, color_type):
    """Calculates the bolometric correction, using a polynomial fit.

    Args:
        color_value (float): B-V, V-I, or B-I color of the supernova in
            magnitudes (corrected for reddening and extinction from the
            host and MWG.)
        color_err (float): Uncertainty in the photometric color.
        color_type (str): String signifying which color color_value represents.
            Valid values are "BminusV" for B-V, "VminusI" for V-I, and
            "BminusI" for B-I.

    Returns:
        tuple: A tuple containing the bolometric correction for use in
        calculating the bolometric luminosity of the supernova, and the
        uncertainty in that bolometric correction (if the color given
        is within the valid range of the polynomial fit.)

        (bolometric_correction, uncertainty)

        (-999, -999) if the color is outside the valid range.
    """
    bolometric_correction = 0.0

    coefficients, range_min, range_max, rms_err = set_constants(color_type)

    if valid_color(color_value, range_min, range_max):
        bolometric_correction = calculate_polynomial(coefficients,
                                                     color_value)
        uncertainty = calc_bolometric_correction_err(color_value,
                                                     color_err, color_type)
    else:
        bolometric_correction = -999
        uncertainty = -999

    return bolometric_correction, uncertainty
