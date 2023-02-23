from superbol.lum import BolometricFlux
from superbol import fqbol
from superbol import planck 
from superbol.blackbody import *
from superbol import mag2flux
from operator import attrgetter

def trim_SED(SED, min_wavelength):
    """
    Trim SED to only include fluxes greater than the `min_wavelength`
    
    Args: 
        SED (list): list of SEDs
        min_wavelength (number): minimum wavelength that fluxes must have
    
    Returns: 
        list: SEDs that meet minimum wavelength

    """
    return [flux for flux in SED if flux.wavelength >= min_wavelength]
    # ^^ Trims everything shorter than the max flux's wavelength
    
# A function to find the max wavelength in the list of SEDs
# 2 readings with SAME FLUX: use shortest wavelength (Max is fine, then. Wavelengths sorted already).
def find_max_flux(SED):
    ''' Find the max flux in list SEDs. '''
    
    max_flux = max(SED, key= attrgetter('flux'))

    # Grabs the max flux's wavelength for return
    # IMPORTANT: Max() returns the FIRST instance it encounters.
    
    return max_flux

def find_min_flux(SED):
    ''' Finds the minimum flux in list SEDs.'''
    min_flux = min(SED, key= attrgetter('flux'))

    return min_flux

def trim_SED_to_peak(SED):
    ''' Trims SED to keep only the fluxes whose wavelengths are greater/equal to that
        of the peak'''

    max_flux = find_max_flux(SED)
    return [flux for flux in SED if flux.wavelength >= max_flux.wavelength]


def get_augmented_bolometric_flux(SED):
    
    # Grabs the max flux from the SED
    max_flux = find_max_flux(SED)

    bbfit = BlackbodyFit() #Blackbody fit object

    bbfit.fit_to_SED(SED)   # Fit the blackbody to the SED 


    min_flux = find_min_flux(SED)
    min_wavelength = min_flux.wavelength
    
    # Calculate the IR correction. (Integration from 0 to ∞) - (0 to max_wavelength)
    # Integrates from max_flux_wavelength to ∞ (Infrared range)

    uv_correction = bb_flux_integrated(min_wavelength, bbfit.temperature, bbfit.angular_radius)

    ir_correction = (bb_total_flux(bbfit.temperature, bbfit.angular_radius) - uv_correction)

    # Integrate straight line from shortest wavelength (original SED, not trimmed) down to 0 at 2k angstroms
    return uv_correction + quasi_bol_flux + ir_correction
