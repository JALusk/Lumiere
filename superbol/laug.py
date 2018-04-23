def trim_SED(SED, min_wavelength=5000):
    """Trim SED to only include fluxes greater than the `min_wavelength`"""
    trimmed_SED = []
    for flux in SED:
        if flux.wavelength > min_wavelength:
            trimmed_SED.append(flux)

    return trimmed_SED
