def trim_SED(SED, min_wavelength=5000):
    """
    Trim SED to only include fluxes greater than the `min_wavelength`
    
    Args: 
        SED (list): list of SEDs
        min_wavelength (number): minimum wavelength that fluxes must have
    
    Returns: 
        list: SEDs that meet minimum wavelength

    """
    return [flux for flux in SED if flux.wavelength > min_wavelength]
