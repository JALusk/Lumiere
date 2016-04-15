import numpy as np
from scipy.optimize import curve_fit
from astropy import units as u
from astropy import constants as const
from planck import *

def bb_flux(wavelength, temperature, angular_radius):
    """Observed flux at [wavelength] from a blackbody of [temperature] and [angular_radius]
    """
    bb_flux = (np.pi * u.sr) * planck_function(wavelength, temperature) * (angular_radius)**2

    return bb_flux

def bb_flux_nounits(wavelength, temperature, angular_radius):
    flux = bb_flux(wavelength, temperature, angular_radius)

    return flux.value

def bb_flux_integrated(wavelength, temperature, angular_radius):
    """Integrate the planck function from 0 to [wavelength], then convert result to observed flux
    """
    bb_flux_integrated = (np.pi * u.sr) * planck_integral(wavelength, temperature) * (angular_radius)**2

    return bb_flux_integrated.value

def dbb_flux_integrated_dT(wavelength, temperature, angular_radius):
    """Take the derivative of the integrated planck function, then convert result to observed flux
    """
    dbb_flux_integrated_dT = (np.pi * u.sr) * d_planck_integral_dT(wavelength, temperature) * (angular_radius)**2

    return dbb_flux_integrated_dT.value

def bb_total_flux(temperature, angular_radius):
    """Integrate the planck function from 0 to infinity, then convert result to observed flux
    """
    temperature = u.Quantity(temperature, u.K)
    bb_total_flux = const.sigma_sb * (angular_radius**2) * temperature**4
    bb_total_flux = bb_total_flux.to(u.erg / (u.s * u.cm**2))
    return bb_total_flux.value

def dbb_total_flux_dT(temperature, angular_radius):
    """Take the derivative of the total blackbody flux with respect to temperature
    """
    temperature = u.Quantity(temperature, u.K)
    bb_total_flux = 4.0 * const.sigma_sb * (angular_radius**2) * temperature**3
    bb_total_flux = bb_total_flux.to(u.erg / (u.s * u.cm**2 * u.K))
    return bb_total_flux.value

def dBB_dT(wavelength, temperature, angular_radius):
    dBB_dT = (np.pi) * dplanck_dT(wavelength, temperature) * (angular_radius)**2

    return dBB_dT

def dBB_dT_nounits(wavelength, temperature, angular_radius):
    dBB_dT_nounits = dBB_dT(wavelength, temperature, angular_radius)

    return dBB_dT_nounits.value

def calculate_chisq(y_data, y_data_uncertainties, x_data, func, parameters):
    chisq = np.sum(((y_data - func(x_data, *parameters))/y_data_uncertainties)**2)
    return chisq
    
def bb_fit_parameters(wavelengths, fluxes, flux_uncertainties):
    popt, pcov = curve_fit(bb_flux_nounits, wavelengths, fluxes, p0=[5000, 1.0e-10])
    temperature = popt[0]
    angular_radius = popt[1]
    perr = np.sqrt(np.diag(pcov))
    
    return temperature, angular_radius, perr
