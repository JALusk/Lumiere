from superbol.lum import BolometricFlux
from superbol import fqbol
from superbol import planck 
from superbol.blackbody import *
from superbol import mag2flux
from operator import attrgetter

def trim_SED(SED, wavelength = 0):
    """
    Trim SED to only include fluxes greater than the `min_wavelength`
    
    Args: 
        SED (list): list of SEDs
        min_wavelength (number): minimum wavelength that fluxes must have
    
    Returns: 
        list: SEDs that meet minimum wavelength

    """
    if not wavelength:
        trimmedSED = trim_SED_to_peak(SED)
    else:
        trimmedSED =  [flux for flux in SED if flux.wavelength >= wavelength]
    return trimmedSED
    
def find_max_flux(SED):
    ''' Find the max flux in list SEDs. '''
    max_flux = max(SED, key= attrgetter('flux'))
    # IMPORTANT: Max() returns the FIRST instance it encounters.
    return max_flux

def find_longest_wavelength(SED):
    list1 = sorted(flux.wavelength for flux in SED)
    return list1[-1]

def find_shortest_wavelength(SED):
    list1 = sorted(flux.wavelength for flux in SED)
    return list1[0]


def find_min_flux(SED):
    ''' Finds the minimum flux in list SEDs.'''
    min_flux = min(SED, key= attrgetter('flux'))

    return min_flux

def trim_SED_to_peak(SED):
    ''' Trims SED to keep only the fluxes whose wavelengths are greater/equal to that
        of the peak'''

    max_flux = find_max_flux(SED)
    return [flux for flux in SED if flux.wavelength >= max_flux.wavelength]

def get_UV(SED):
    shortest_wavelength = find_shortest_wavelength(SED)
    bluest_flux = find_min_flux(SED)
    
    uv_correction = 0.5 * (shortest_wavelength * bluest_flux.flux)
    return uv_correction

def get_IR(SED):
    longest_wavelength = find_longest_wavelength(SED)

    bbfit = BlackbodyFit()
    trimmedSED = trim_SED(SED)
    bbfit.fit_to_SED(trimmedSED) 

    f_ir_trimmed = (bb_total_flux(bbfit.temperature, bbfit.angular_radius) 
                    - bb_flux_integrated(longest_wavelength, bbfit.temperature, bbfit.angular_radius))

    return f_ir_trimmed

def get_augmented_bolometric_flux(SED):
    pass

