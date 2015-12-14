import numpy as np
from scipy.optimize import curve_fit
from astropy import units as u
from planck import planck_function, dplanck_dT

def bb_flux(wavelength, temperature, angular_radius):
    bb_flux = (np.pi * u.sr) * planck_function(wavelength, temperature) * (angular_radius)**2

    return bb_flux

def bb_flux_nounits(wavelength, temperature, angular_radius):
    flux = bb_flux(wavelength, temperature, angular_radius)

    return flux.value

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
