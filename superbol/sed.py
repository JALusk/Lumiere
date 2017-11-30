import math
import numpy as np

from itertools import groupby
from superbol import mag2flux

def group_fluxes(fluxes, keyfunc=math.floor):
    """Group fluxes by applying keyfunc() to each flux.time in the fluxes.

    Fluxes with the same return value for keyfunc(flux.time) will be 
    part of the same group."""
    grouped_fluxes = [list(it) for k, it in groupby(fluxes, lambda x: keyfunc(x.time))]
    return grouped_fluxes

def get_weights(uncertainties):
    """Calculate weights from uncertainties"""
    weights = [1/s**2 for s in uncertainties]
    return weights

def weighted_average(values, uncertainties):
    """Calculate the weighed average of values with uncertainties"""
    weights = get_weights(uncertainties)
    denominator = sum(weights)
    numerator = sum([weights[i] * values[i] for i in range(len(values))])
    return numerator/denominator

def weighted_average_uncertainty(uncertainties):
    """Calculate the uncertainty in a weighted average"""
    weights = get_weights(uncertainties)
    uncertainty = 1.0/math.sqrt(sum(weights))
    return uncertainty

def get_flux_values(fluxes):
    return [f.flux for f in fluxes]

def get_flux_uncertainties(fluxes):
    return [f.flux_uncertainty for f in fluxes]

def combine_fluxes(fluxes):
    """Combine a list of MonochromaticFluxes using a weighted average"""
    values = get_flux_values(fluxes)
    uncertainties = get_flux_uncertainties(fluxes)
    combined_flux = weighted_average(values, uncertainties)
    combined_uncertainty = weighted_average_uncertainty(uncertainties)
    wavelength = fluxes[0].wavelength
    time = fluxes[0].time
    return mag2flux.MonochromaticFlux(combined_flux,
                                      combined_uncertainty,
                                      wavelength,
                                      time)

def yield_fluxes_at_each_observed_wavelength(fluxes):
    """Yield lists of MonochromaticFluxes with the same wavelength"""
    for wavelength in set(f.wavelength for f in fluxes):
        yield [f for f in fluxes if f.wavelength == wavelength]

def get_integrable_fluxes(fluxes):
    """Return a list of MonochromaticFluxes with duplicates averaged"""
    integrable_fluxes = []
    for f in yield_fluxes_at_each_observed_wavelength(fluxes):
        if len(f) == 1:
            integrable_fluxes.append(f[0])
        else:
            integrable_fluxes.append(combine_fluxes(f))
    return integrable_fluxes
