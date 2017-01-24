import numpy as np
from astropy import units as u


def mag2flux(magnitude, uncertainty, effective_wl, flux_at_zero_mag):
    """Converts an observed magnitude in a filter band to an average flux.

    Args:
        magnitude (float):  Apparent magnitude.
        uncertainty (float): Apparent magnitude uncertainty.
        effective_wl (float): Effective wavelength of the filter.
        flux_at_zero_mag (float): Flux at zero magnitude of the filter.

    Returns:
        tuple: A tuple of two floats:

        * the flux in :math:`erg \\; s^{-1} cm^{-2} Angstrom^{-1}`
        * the flux uncertainty in :math:`erg \\; s^{-1} cm^{-2} Angstrom^{-1}`

        (flux, flux_uncertainty)
    """
    effective_wl = effective_wl * u.AA
    flux_at_zero_mag = flux_at_zero_mag * (u.erg / (u.s * u.cm**2 * u.AA))

    flux = flux_at_zero_mag * 10**(-0.4 * magnitude)
    flux_uncertainty = np.abs(flux * -0.4 * np.log(10) * uncertainty)

    return flux.value, flux_uncertainty.value
