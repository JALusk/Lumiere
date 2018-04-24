import math
import numpy as np

from itertools import groupby
from superbol import mag2flux

def group_magnitudes(magnitudes, keyfunc=math.floor):
    """Group magnitudes by applying keyfunc() to each magnitude.time.

    Magnitudes with the same return value for keyfunc(magnitude.time) will be 
    part of the same group."""
    grouped_magnitudes = [list(it) for k, it in groupby(magnitudes, lambda x: keyfunc(x.time))]
    return grouped_magnitudes

def combine_magnitudes(observed_magnitudes):
    """Combine a list of ObservedMagnitudes using a weighted average"""
    values = get_observed_magnitude_values(observed_magnitudes)
    uncertainties = get_observed_magnitude_uncertainties(observed_magnitudes)
    combined_magnitude = weighted_average(values, uncertainties)
    combined_uncertainty = weighted_average_uncertainty(uncertainties)
    band = observed_magnitudes[0].band
    time = observed_magnitudes[0].time
    return mag2flux.ObservedMagnitude(combined_magnitude,
                                      combined_uncertainty,
                                      band,
                                      time)

def get_observed_magnitude_values(observed_magnitudes):
    return [obs.magnitude for obs in observed_magnitudes]

def get_observed_magnitude_uncertainties(observed_magnitudes):
    return [obs.uncertainty for obs in observed_magnitudes]

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

def get_weights(uncertainties):
    """Calculate weights from uncertainties"""
    weights = [1/s**2 for s in uncertainties]
    return weights

