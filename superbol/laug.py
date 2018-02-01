def trim_sed(SED, min_wavelength):
    """Trim SED to only include fluxes greater than the `min_wavelength`"""
    trimmed_sed = []
    for flux in SED:
        if flux.wavelength > min_wavelength:
            trimmed_sed.append(flux)

    return trimmed_sed
