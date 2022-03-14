import numpy as np
from astropy import constants as const
from astropy import units as u
import functools


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
    iterations = min(int(2.0 + 20.0/x.value), 512)

    series = functools.reduce(
        lambda acc, i: acc + (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6.0/i**4) * np.exp(-i * x), range(1, iterations), 0)

    B_integral = (C1 * temperature**4 / C2**4) * series / u.sr

    return B_integral.to(u.erg / (u.s * u.cm**2 * u.sr))
