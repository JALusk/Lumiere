import numpy as np
from astropy import units as u

def get_filter_parameters(filter_band):
    """Fetches effective wavelength and flux zeropoint of a filter.

    Args:
        filter_band: StringType - Name of the filter band (ex: 'U').

    Returns:
        A tuple containing two astropy quantities: the effective wavelength
        in Angstroms, and the flux at zero magnitudes in erg/s/cm^s/A.

        (effective_wl, flux_at_zero_mag).

    Raises:
        TypeError: The argument given is not a string.
        KeyError: The argument given is not a valid filter name.
    """
    if not isinstance(filter_band, str):
        raise TypeError('filter_band must be StringType')
    elif filter_band in filter_data:
        effective_wl = filter_data[filter_band]['effective_wl']
        flux_at_zero_mag = filter_data[filter_band]['flux_at_zero_mag']
    else:
        raise KeyError('Your specified filter (%s) is not available', filter_band)
   
    # Attach CGS units to parameters
    effective_wl = effective_wl * u.AA
    flux_at_zero_mag = flux_at_zero_mag * (u.erg / (u.s * u.cm**2 * u.AA)) 

    return effective_wl, flux_at_zero_mag

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
