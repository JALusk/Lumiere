import math

from superbol import sed

def calculate_lightcurve(fluxes, distance, luminosity_calculator):
    grouped_fluxes = sed.group_fluxes(fluxes, math.floor)
    lightcurve = []
    for flux_group in grouped_fluxes:
        if len(flux_group) > 2:
            luminosity = luminosity_calculator(flux_group, distance)
            lightcurve.append(luminosity)
    return lightcurve
