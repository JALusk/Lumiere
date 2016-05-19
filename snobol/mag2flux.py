import numpy as np
from astropy import units as u

def mag2flux(magnitude, uncertainty, effective_wl, flux_at_zero_mag):
    """Converts an observed magnitude in a filter band to an average flux.

    Args:
        magnitude: FloatType - Apparent magnitude.
        uncertainty: FloatType - Apparent magnitude uncertainty.
        effective_wl: FloatType - Effective wavelength of the filter.
        flux_at_zero_mag: FloatType - Flux at zero magnitude of the filter.

    Returns:
        a tuple of two astropy quantities:
            * the flux in erg/s/cm^2/A, 
            * the flux uncertainty in erg/s/cm^2/A

        (flux, flux_uncertainty, effective_wl)
    """
    effective_wl = effective_wl * u.AA
    flux_at_zero_mag = flux_at_zero_mag * (u.erg / (u.s * u.cm**2 * u.AA)) 
   
    flux = flux_at_zero_mag * 10**(-0.4 * magnitude)
    flux_uncertainty = np.abs(flux * -0.4 * np.log(10) * uncertainty)

    return flux.value, flux_uncertainty.value
