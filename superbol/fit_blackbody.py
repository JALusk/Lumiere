import numpy as np
from scipy.optimize import curve_fit
from astropy import units as u
from astropy import constants as const
from superbol.planck import *

def bb_flux(wavelength, temperature, angular_radius):
    """Observed flux at `wavelength` from a blackbody of `temperature` and `angular_radius` in cgs units

    Args:
        wavelength (float): Wavelength in Angstroms
        temperature (float): Temperature in Kelvin
        angular_radius (float): Angular radius :math:`(\\theta = \\frac{R}{D})`

    Returns:
        Astropy Quantity: blackbody flux in :math:`erg \\; s^{-1} cm^{-2} Anstrom^{-1}`
    """
    bb_flux = (np.pi * u.sr) * planck_function(wavelength, temperature) * (angular_radius)**2

    return bb_flux

def bb_flux_nounits(wavelength, temperature, angular_radius):
    """Identical to bb_flux, but returns the value rather than the Astropy Quantity.

    This function is needed because curve_fit() does not play nice with functions that return Astropy quantities.
    """
    flux = bb_flux(wavelength, temperature, angular_radius)

    return flux.value

def bb_flux_integrated(wavelength, temperature, angular_radius):
    """Integrate the planck function from :math:`\\lambda = 0` to :math:`\\lambda =` `wavelength`, then convert result to observed flux

    Args:
        wavelength (float): Wavelength in Angstroms
        temperature (float): Temperature in Kelvin
        angular_radius (float): Angular radius :math:`(\\theta = \\frac{R}{D})`

    Returns:
        float: Value of the integrated flux in :math:`erg \\; s^{-1} cm^{-2}`
    """
    bb_flux_integrated = (np.pi * u.sr) * planck_integral(wavelength, temperature) * (angular_radius)**2

    return bb_flux_integrated.value

def dbb_flux_integrated_dT(wavelength, temperature, angular_radius):
    """Take the derivative of the integrated planck function, then convert result to observed flux. This is used in error propagation calculations.

    Args:
        wavelength (float): Wavelength in Angstroms
        temperature (float): Temperature in Kelvin
        angular_radius (float): Angular radius :math:`(\\theta = \\frac{R}{D})`

    Returns:
        float: Derivative of the integrated blackbody flux at `wavelength` with respect to `temperature`
    """
    dbb_flux_integrated_dT = (np.pi * u.sr) * d_planck_integral_dT(wavelength, temperature) * (angular_radius)**2

    return dbb_flux_integrated_dT.value

def bb_total_flux(temperature, angular_radius):
    """Integrate the planck function from :math:`\\lambda = 0` to :math:`\\lambda = \\infty`, then convert result to observed flux.

    The calculation is done using the Stefan-Boltxmann law, rather than performing an actual integral using numerical integration. The `angular_radius` is included to convert the result to an observed flux.

    Args:
        temperature (float): Temperature in Kelvin
        angular_radius (float): Angular radius :math:`(\\theta = \\frac{R}{D})`

    Returns:
        float: The value of :math:`\\sigma \\frac{R^2}{D^2} T^4`
    """
    temperature = u.Quantity(temperature, u.K)
    bb_total_flux = const.sigma_sb * (angular_radius**2) * temperature**4
    bb_total_flux = bb_total_flux.to(u.erg / (u.s * u.cm**2))
    return bb_total_flux.value

def dbb_total_flux_dT(temperature, angular_radius):
    """Take the derivative of the total blackbody flux with respect to temperature

    The calculation takes a derivative of the Stefan-Boltzmann law with respect to temperature. The `angular_radius` is present to convert the result to an observed flux.

    Args:
        temperature (float): Temperature in Kelvin
        angular_radius (float): Angular radius :math:`(\\theta = \\frac{R}{D})`

    Returns:
        float: The value of :math:`4 \\sigma \\frac{R^2}{D^2} T^3`
    """
    temperature = u.Quantity(temperature, u.K)
    bb_total_flux = 4.0 * const.sigma_sb * (angular_radius**2) * temperature**3
    bb_total_flux = bb_total_flux.to(u.erg / (u.s * u.cm**2 * u.K))
    return bb_total_flux.value

def bb_fit_parameters(wavelengths, fluxes, flux_uncertainties):
    """Fit a blackbody to observed `wavelengths` and `fluxes`.

    The initial guesses for the `temperature` and `angular_radius` are :math:`T = 5000` K and :math:`\\theta = 1.0 \\times 10^{-10}`. These are typical for an extragalactic supernovae, but should be used with caution for any other objects.

    Args:
        wavelengths (list): List of wavelengths at which the fluxes were observed.
        fluxes (list): List of observed monochromatic fluxes.

    Returns:
        tuple: Tuple containing the best fit `temperature`, `angular_radius`, as well as perr, a 2-tuple containing the `temperature` error and the `angular_radius` error.

        (temperature, angular_radius, perr)
    """
    popt, pcov = curve_fit(bb_flux_nounits, wavelengths, fluxes, p0=[5000, 1.0e-10], sigma=flux_uncertainties, absolute_sigma=True)
    temperature = popt[0]
    angular_radius = popt[1]
    perr = np.sqrt(np.diag(pcov))

    return temperature, angular_radius, perr
