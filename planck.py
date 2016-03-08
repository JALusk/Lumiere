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
