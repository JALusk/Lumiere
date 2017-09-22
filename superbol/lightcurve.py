import math

from itertools import groupby
from superbol import lqbol


def group_fluxes(fluxes, keyfunc=None):
    """Group fluxes by applying keyfunc() to each flux.time in the fluxes.

    Fluxes with the same return value for keyfunc(flux.time) will be 
    part of the same group."""
    grouped_fluxes = [list(it) for k, it in groupby(fluxes, lambda x: keyfunc(x.time))]
    return grouped_fluxes

def calculate_lightcurve(fluxes, distance, luminosity_calculator):
    grouped_fluxes = group_fluxes(fluxes, math.floor)
    lightcurve = []
    for flux_group in grouped_fluxes:
        if len(flux_group) > 2:
            luminosity = luminosity_calculator(flux_group, distance)
            lightcurve.append(luminosity)
    return lightcurve
