import numpy as np
from astropy import units as u
import scipy.integrate as integrate
from mag2flux import *
from specutils import extinction
from fit_blackbody import *

def integrate_fqbol(wavelengths, fluxes, flux_uncertainties):
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
    ir_correction = integrate.quad(bb_flux_nounits, longest_wl, np.inf,
                                   args=(temperature, angular_radius))[0]

    T_errterm = integrate.quad(dBB_dT_nounits, longest_wl, np.inf,
                               args=(temperature, angular_radius))[0] * T_err
    rad_errterm = 2 * ir_correction / angular_radius * rad_err

    ir_corr_err = np.sqrt(T_errterm**2 + rad_errterm**2)

    return ir_correction, ir_corr_err

def uv_correction_blackbody(temperature, T_err, angular_radius, rad_err, shortest_wl):
    uv_correction = integrate.quad(bb_flux_nounits, 0.0, shortest_wl, 
                                   args=(temperature, angular_radius))[0]

    T_errterm = integrate.quad(dBB_dT_nounits, 500.0, shortest_wl,
                               args=(temperature, angular_radius))[0] * T_err
    rad_errterm = 2 * uv_correction / angular_radius * rad_err

    uv_corr_err = np.sqrt(T_errterm**2 + rad_errterm**2)

    return uv_correction, uv_corr_err

def uv_correction_linear(shortest_wl, shortest_flux, shortest_flux_err):
    fluxes = [0.0, shortest_flux]
    wavelengths = [2000.0, shortest_wl]
    uv_correction = np.trapz(fluxes, wavelengths)
    uv_correction_err = 0.5 * (shortest_wl - 2000.0) * shortest_flux_err
    
    return uv_correction, uv_correction_err
