import math

from superbol import sed

def calculate_lightcurve(fluxes, distance, luminosity_calculator):
    grouped_fluxes = sed.group_fluxes(fluxes, math.floor)
    lightcurve = []
    SEDs = []
    for flux_group in grouped_fluxes:
        if len(flux_group) > 2:
            SEDs.append(sed.get_SED(flux_group))

    sed.interpolate_missing_fluxes(SEDs)

    for SED in SEDs:
        luminosity = luminosity_calculator(SED, distance)
        lightcurve.append(luminosity)
    return lightcurve
