import numpy as np
from astropy import constants as const
from astropy import units as u

def planck_function(wavelength, temperature):
    """Planck function at given wavelength and temperature in cgs.

    Args:
        wavelength: Wavelength in Angstrom
        tempereature: Temperature in Kelvin
    
    Returns:
        B_lambda: The specific intensity of the Planck function in
                  erg s^-1 cm^-2 sterad^-1 AA^-1
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
    """ Integrate the Planck function from lambda = 0 to lambda = wavelength
        using the infinite series approximation of the integral.

    Args:
        wavelength: Upper bound for the wavelength in Angstrom.
        temperature: Temperature in Kelvin.
    Returns:
        B_integral: Integral of the planck function
    """
    wavelength = u.Quantity(wavelength, unit=u.Angstrom)
    temperature = u.Quantity(temperature, unit=u.K)

    C1 = 2.0 * const.h.cgs * const.c.cgs**2
    C2 = const.h.cgs * const.c.cgs / const.k_B.cgs

    x = C2 / (wavelength.to(u.cm) * temperature)
    iterations = min(int(2 + 20/x.value), 512)

    series = 0.0
    for i in range(1, iterations):
        term = (x**3/i + 3*x**2/i**2 + 6*x/i**3 + 6/i**4) * np.exp(-i * x)
        series += term

    B_integral = (C1 * temperature**4 / C2**4) * series / u.sr
    return B_integral.to(u.erg / (u.s * u.cm**2 * u.sr))

def d_planck_integral_dT(wavelength, temperature):
    """ Derivative of the integrated Planck function from lambda = 0 to lambda = wavelength using the infinite series approximation of the integral. This is used in the error propagation

    Args:
        wavelength: Upper bound for the wavelength in Angstrom.
        temperature: Temperature in Kelvin.
    Returns:
        dB_integral_dT: Derivative of the integral of the planck function with respect to T
    """
    wavelength = u.Quantity(wavelength, unit=u.Angstrom)
    temperature = u.Quantity(temperature, unit=u.K)

    C1 = 2.0 * const.h.cgs * const.c.cgs**2
    C2 = const.h.cgs * const.c.cgs / const.k_B.cgs

    x = C2 / (wavelength.to(u.cm) * temperature)
    iterations = min(int(2 + 20/x.value), 512)

    series = 0.0
    for i in range(1, iterations):
        term1 = C1 / (temperature * wavelength**4)
        term2 = 4*C1 / (i * C2 * wavelength**3)
        term3 = 12 * C1 * temperature / (i**2 * C2**2 * wavelength**2)
        term4 = 24*C1*temperature**2 / (i**3 * C2**3 * wavelength)
        term5 = 24 * C1 * temperature**3 / (i**4 * C2**4)
        series += (term1 + term2 + term3 + term4 + term5) * np.exp(-i*x)

    dB_integral_dT = series / u.sr
    return dB_integral_dT.to(u.erg / (u.s * u.cm**2 * u.K * u.sr))

def dplanck_dT(wavelength, temperature):
    """Partial derivative of the Planck function with Temperature

    Args:
        wavelength: Wavelength in Angstrom
        temperature: Temperature in Kelvin

    Returns:
        dB_lambda_dT: derivative of the specific intensity of the Planck function in erg s^-1 cm^-2 sterad^-1 AA^-1 K^-1
    """
    wavelength = u.Quantity(wavelength, unit=u.Angstrom)
    temperature = u.Quantity(temperature, unit=u.K)

    C1 = 2.0 * const.h.cgs * const.c.cgs**2
    C2 = const.h.cgs * const.c.cgs / const.k_B.cgs

    exp_term = np.exp(C2 / (wavelength * temperature))
    expm1_term = np.expm1(C2 / (wavelength * temperature))

    dB_lambda_dT = (C1 * C2 * exp_term) / (wavelength**6 * temperature**2 * (expm1_term)**2) / u.sr
    dB_lambda_dT = dB_lambda_dT.to(u.erg / (u.s * u.cm**2 * u.AA * u.sr * u.K))

    return dB_lambda_dT
