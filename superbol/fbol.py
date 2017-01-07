import numpy as np

from superbol.fit_blackbody import (bb_flux_integrated, bb_total_flux,
                                    dbb_flux_integrated_dT, dbb_total_flux_dT)


def integrate_fqbol(wavelengths, fluxes, flux_uncertainties):
    """Calculate the trapezoidal rule integral of the observed `fluxes`.

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
            term = 0.5 * (wavelengths[i + 1] - wavelengths[i]) * uncertainty
            quad_terms = np.append(quad_terms, term)
        elif i == len(flux_uncertainties) - 1:
            term = 0.5 * (wavelengths[i] - wavelengths[i - 1]) * uncertainty
            quad_terms = np.append(quad_terms, term)
        else:
            term = 0.5 * (
                wavelengths[i + 1] - wavelengths[i - 1]) * uncertainty
            quad_terms = np.append(quad_terms, term)
    fqbol_uncertainty = np.sqrt(np.sum(x * x for x in quad_terms))

    fqbol_uncertainty = fqbol_uncertainty

    return fqbol, fqbol_uncertainty


def ir_correction(temperature, T_err, angular_radius, rad_err, longest_wl):
    """Apply correction for unobserved flux in the IR.

    After the temperature and angular radius has been found through fitting a
    blackbody to the observed fluxes, this function takes those values and
    integrates under the fitted blackbody function from the longest observed
    wavelength out to :math:`\\lambda = \\infty`.

    Args:
        temperature (float): Best fit blackbody temperature in Kelvin
        T_err (float): Uncertainty in best fit blackbody temperature in Kelvin
        angular_radius (float): Best fit blackbody angular radius
        rad_err (float): Uncertainty in best fit blackbody angular radius
        longest_wl (float): Longest observed wavelength

    Returns:
        tuple: 2-tuple

        * (float): The IR correction in :math:`erg \\; s^{-1} cm^{-2}`
        * (float): The uncertainty in the IR correction in the same units
    """
    ir_correction = bb_total_flux(temperature,
                                  angular_radius) - bb_flux_integrated(
                                      longest_wl, temperature, angular_radius)

    T_errterm = (dbb_total_flux_dT(temperature, angular_radius) -
                 dbb_flux_integrated_dT(longest_wl, temperature,
                                        angular_radius)) * T_err
    rad_errterm = 2 * ir_correction / angular_radius * rad_err

    ir_corr_err = np.sqrt(T_errterm**2 + rad_errterm**2)

    return ir_correction, ir_corr_err


def uv_correction_blackbody(temperature, T_err, angular_radius, rad_err,
                            shortest_wl):
    """Apply correction for unobserved flux in the UV using the blackbody fit.

    After the temperature and angular radius have been found through fitting a
    blackbody to the observed fluxes, this function takes those values and
    integrates under the fitted blackbody from the shortest observed wavelength
    down to :math:`\\lambda = 0`.

    Args:
        temperature (float): Best fit blackbody temperature in Kelvin
        T_err (float): Uncertainty in best fit blackbody temperature in Kelvin
        angular_radius (float): Best fit blackbody angular radius
        rad_err (float): Uncertainty in best fit blackbody angular radius
        shortest_wl (float): Shortest observed wavelength

    Returns:
        tuple: 2-tuple

        * (float): The UV correction in :math:`erg \\; s^{-1} cm^{-2}`
        * (float): The uncertainty in the UV correction in the same units
    """
    uv_correction = bb_flux_integrated(shortest_wl, temperature,
                                       angular_radius)

    T_errterm = dbb_flux_integrated_dT(shortest_wl, temperature,
                                       angular_radius) * T_err
    rad_errterm = 2 * uv_correction / angular_radius * rad_err

    uv_corr_err = np.sqrt(T_errterm**2 + rad_errterm**2)

    return uv_correction, uv_corr_err


def uv_correction_linear(shortest_wl, shortest_flux, shortest_flux_err):
    """Apply correction for unobserved flux in the UV using a linear function.

    This function integrates under a straight line from the shortest observed
    wavelength down to :math:`f(\\lambda) = 0` at :math:`\\lambda = 2000`
    Angstroms. This approximates the effects of line blanketing in the UV as in
    Bersten & Hamuy (2009).

    Args:
        shortest_wl (float): Shortest observed wavelength
        shortest_flux (float): Flux at shortest observed wavelength
        shortest_flux_err (float): Uncertainty in the shortest observed flux

    Returns:
        tuple: 2-tuple

        * (float): The UV correction in :math:`erg \\; s^{-1} cm^{-2}`
        * (float): The uncertainty in the UV correction in the same units
    """
    fluxes = [0.0, shortest_flux]
    wavelengths = [2000.0, shortest_wl]
    uv_correction = np.trapz(fluxes, wavelengths)
    uv_correction_err = 0.5 * (shortest_wl - 2000.0) * shortest_flux_err

    return uv_correction, uv_correction_err
