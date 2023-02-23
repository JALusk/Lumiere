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
    trimmedSED =  [flux for flux in SED if flux.wavelength >= min_wavelength]
    return trimmedSED
    # ^^ Trims everything shorter than the max flux's wavelength
    
# A function to find the max wavelength in the list of SEDs
# 2 readings with SAME FLUX: use shortest wavelength (Max is fine, then. Wavelengths sorted already).
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
    # UV correction = 1/2(B*H)
    max_flux = find_max_flux(SED)
    longest_wavelength = max_flux.wavelength
    base = int(longest_wavelength) - 2000
    height = find_min_flux(SED)
    uv_correction = (0.5 * (base * int(height.flux)))

    return uv_correction

def get_IR(SED):
    longest_wavelength = find_longest_wavelength(SED)
    bbfit = BlackbodyFit()  # Blackbody object
    trimmedSED = trim_SED(SED, 5000)
    bbfit.fit_to_SED(trimmedSED)   # fits the blackbody to the SED

    f_ir_trimmed = (bb_total_flux(bbfit.temperature, bbfit.angular_radius) 
                    - bb_flux_integrated(longest_wavelength, bbfit.temperature, bbfit.angular_radius))

    return f_ir_trimmed





def get_augmented_bolometric_flux(SED):


    bbfit = BlackbodyFit() #Blackbody fit object

    trimmedSED = trim_SED(SED, 5000)
    bbfit.fit_to_SED(trimmedSED)   # Fit the blackbody to the SED 


    min_flux = find_min_flux(SED)
    min_wavelength = min_flux.wavelength

    # Calculate the IR correction. (Integration from 0 to ∞) - (0 to max_wavelength)
    # Integrates from max_flux_wavelength to ∞ (Infrared range)

    #hardly ever need to integrate UV. If statement (rarely happens)
    uv_correction = bb_flux_integrated(min_wavelength, bbfit.temperature, bbfit.angular_radius)

    #UV: base: longest wavelength - 2000, height is flux (1/2BH)
    # Get 0 to infinity, subtract before longest wavelength
    ir_correction = (bb_total_flux(bbfit.temperature, bbfit.angular_radius) - uv_correction)

    # Integrate straight line from shortest wavelength (original SED, not trimmed) down to 0 at 2k angstroms


    # To get quasi-bol-flux:
    #   Call out to quasi bolometric function and get result

    #methods in faug returns IR number and UV number (IR correction function, UV correction function)
    # dont worry about augmented bol flux yet. That'll go in a different file/hierarchy

    return uv_correction + quasi_bol_flux + ir_correction
