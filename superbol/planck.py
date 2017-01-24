import numpy as np
from astropy import constants as const
from astropy import units as u


def planck_function(wavelength, temperature):
    """Planck function at given `wavelength` and `temperature` in cgs.

    The specific intensity of radiation from a blackbody source is given by:

    :math:`B_{\\lambda}(\\lambda, T) = \\displaystyle\\frac{C_1}{\\lambda^5} \\frac{1}{e^{\\frac{C_2}{\\lambda T}} - 1}`

    where :math:`C_1 = 2hc^2` and :math:`C_2 = hc/k_B` are the first and second radiation constants.

    Args:
        wavelength (float): Wavelength in Angstrom
        tempereature (float): Temperature in Kelvin

    Returns:
        Astropy Quantity: The specific intensity of the
        Planck function in :math:`erg \\; s^{-1} cm^{-2} sterad^{-1} Angstrom^{-1}`
    """
    wavelength = u.Quantity(wavelength, unit=u.Angstrom)
    temperature = u.Quantity(temperature, unit=u.K)

    C1 = 2.0 * const.h.cgs * const.c.cgs**2
    C2 = const.h.cgs * const.c.cgs / const.k_B.cgs

    B_lambda = (C1 / wavelength**5) / \
               (np.expm1(C2 / (wavelength * temperature))) / u.sr
    B_lambda = B_lambda.to(u.erg / (u.s * u.cm**2 * u.AA * u.sr))

    return B_lambda


def planck_integral(wavelength, temperature):
    """Integrate the Planck function over a finite wavelength interval.

    The integral is taken from :math:`\\lambda = 0` to :math:`\\lambda =` `wavelength`
    using the infinite series approximation of the integral:

    :math:`\\displaystyle\\int_{0}^{\\lambda_1} B_{\\lambda}(\\lambda, T) \\; d\\lambda = \\frac{C_1 T^4}{C_2^4} \\sum_{n = 1}^{\\infty}\\left(\\frac{x_1^3}{n} + \\frac{3x_1^2}{n^2} + \\frac{6x_1}{n^3} + \\frac{6}{n^4}\\right) e^{-nx_1}`

    where :math:`C_1 = 2hc^2` and :math:`C_2 = hc/k_B` are the first and second radiation constants, and :math:`x_1 = C_2/\\lambda_1 T`
    is a change of variables to make the integral easier to represent in this series form.

    Args:
        wavelength (float): Upper bound for the wavelength in Angstrom.
        temperature (float): Temperature in Kelvin.

    Returns:
        Astropy Quantity: Integral of the planck function in
        :math:`erg \\; s^{-1} cm^{-2} sterad^{-1}`
    """
    wavelength = u.Quantity(wavelength, unit=u.Angstrom)
    temperature = u.Quantity(temperature, unit=u.K)

    C1 = 2.0 * const.h.cgs * const.c.cgs**2
    C2 = const.h.cgs * const.c.cgs / const.k_B.cgs

    x = C2 / (wavelength.to(u.cm) * temperature)
    iterations = min(int(2 + 20 / x.value), 512)

    series = 0.0
    for i in range(1, iterations):
        term = (x**3 / i + 3 * x**2 / i**2 + 6 * x / i**3 + 6 / i**4) * np.exp(
            -i * x)
        series += term

    B_integral = (C1 * temperature**4 / C2**4) * series / u.sr

    return B_integral.to(u.erg / (u.s * u.cm**2 * u.sr))


def d_planck_integral_dT(wavelength, temperature):
    """Derivative of the integrated Planck function from :math:`\\lambda = 0` to
    :math:`\\lambda =` `wavelength` using the infinite series approximation of the
    integral. This is used in the error propagation calculation.

    Args:
        wavelength (float): Upper bound for the wavelength in Angstrom.
        temperature (float): Temperature in Kelvin.

    Returns:
        Astropy Quantity: Derivative of the integral of the planck function
        with respect to T in :math:`erg \\; s^{-1} cm^{-2} sterad^{-1} K^{-1}`
    """
    wavelength = u.Quantity(wavelength, unit=u.Angstrom)
    temperature = u.Quantity(temperature, unit=u.K)

    C1 = 2.0 * const.h.cgs * const.c.cgs**2
    C2 = const.h.cgs * const.c.cgs / const.k_B.cgs

    x = C2 / (wavelength.to(u.cm) * temperature)
    iterations = min(int(2 + 20 / x.value), 512)

    series = 0.0
    for i in range(1, iterations):
        term1 = C1 / (temperature * wavelength**4)
        term2 = 4 * C1 / (i * C2 * wavelength**3)
        term3 = 12 * C1 * temperature / (i**2 * C2**2 * wavelength**2)
        term4 = 24 * C1 * temperature**2 / (i**3 * C2**3 * wavelength)
        term5 = 24 * C1 * temperature**3 / (i**4 * C2**4)
        series += (term1 + term2 + term3 + term4 + term5) * np.exp(-i * x)

    dB_integral_dT = series / u.sr

    return dB_integral_dT.to(u.erg / (u.s * u.cm**2 * u.K * u.sr))
