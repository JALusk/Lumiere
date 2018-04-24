import math

from itertools import groupby

def group_magnitudes(fluxes, keyfunc=math.floor):
    """Group fluxes by applying keyfunc() to each flux.time in the fluxes.

    Fluxes with the same return value for keyfunc(flux.time) will be 
    part of the same group."""
    grouped_fluxes = [list(it) for k, it in groupby(fluxes, lambda x: keyfunc(x.time))]
    return grouped_fluxes
