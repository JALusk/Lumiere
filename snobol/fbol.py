import numpy as np
from astropy import units as u
import scipy.integrate as integrate
from mag2flux import *
from specutils import extinction
from fit_blackbody import *

def integrate_fqbol(wavelengths, fluxes, flux_uncertainties):
    """Calculate the trapezoidal rule integral of the observed `fluxes`, and the uncertainty in that integration.
    
    The trapezoidal rule integrates the data by assuming the function is linear between observed points, and then integrates under those line segments.
    The numpy function `trapz` is used to perform the integration, but the uncertainty in the integral due to uncertainties in the observed flux is calculated by hand using standard error propagation techniques.

    Args:
        wavelengths (list): List of wavelengths at which the flux was observed.
        fluxes (list): List of observed fluxes.
        flux_uncertainties (list): List of uncertainties in each observed flux.

    Returns:
        tuple: 2-tuple of floats.

        * The value of the integral
        * The uncertainty in the integral due to uncertainties in the fluxes.

        (fqbol, fqbol_uncertainty)
    """
    fqbol = np.trapz(fluxes, wavelengths)

    quad_terms = np.array([])

    for i, uncertainty in enumerate(flux_uncertainties):
        if i == 0:
            term = 0.5 * (wavelengths[i+1] - wavelengths[i]) * uncertainty
            quad_terms = np.append(quad_terms, term)
        elif i == len(flux_uncertainties) - 1:
            term = 0.5 * (wavelengths[i] - wavelengths[i-1]) * uncertainty
            quad_terms = np.append(quad_terms, term)
        else:
            term = 0.5 * (wavelengths[i+1] - wavelengths[i-1]) * uncertainty
            quad_terms = np.append(quad_terms, term)
    fqbol_uncertainty = np.sqrt(np.sum(x*x for x in quad_terms))

    fqbol_uncertainty = fqbol_uncertainty

    return fqbol, fqbol_uncertainty

def ir_correction(temperature, T_err, angular_radius, rad_err, longest_wl):
    """
    """
    ir_correction = bb_total_flux(temperature, angular_radius) - bb_flux_integrated(longest_wl, temperature, angular_radius)

    T_errterm = (dbb_total_flux_dT(temperature, angular_radius) - dbb_flux_integrated_dT(longest_wl, temperature, angular_radius)) * T_err
    rad_errterm = 2 * ir_correction / angular_radius * rad_err

    ir_corr_err = np.sqrt(T_errterm**2 + rad_errterm**2)

    return ir_correction, ir_corr_err

def uv_correction_blackbody(temperature, T_err, angular_radius, rad_err, shortest_wl):
    uv_correction = bb_flux_integrated(shortest_wl, temperature, angular_radius)

    T_errterm = dbb_flux_integrated_dT(shortest_wl, temperature, angular_radius)* T_err
    rad_errterm = 2 * uv_correction / angular_radius * rad_err

    uv_corr_err = np.sqrt(T_errterm**2 + rad_errterm**2)

    return uv_correction, uv_corr_err

def uv_correction_linear(shortest_wl, shortest_flux, shortest_flux_err):
    fluxes = [0.0, shortest_flux]
    wavelengths = [2000.0, shortest_wl]
    uv_correction = np.trapz(fluxes, wavelengths)
    uv_correction_err = 0.5 * (shortest_wl - 2000.0) * shortest_flux_err
    
    return uv_correction, uv_correction_err
