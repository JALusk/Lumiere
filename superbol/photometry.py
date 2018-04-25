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

def combine_observed_magnitudes(observed_magnitudes):
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

def yield_observed_magnitudes_at_each_observed_band(observed_magnitudes):
    """Yield lists of ObservedMagnitudes with the same band"""
    list_of_bands = sorted(list(set(obs.band for obs in observed_magnitudes)))
    for band in list_of_bands:
        yield [obs for obs in observed_magnitudes if obs.band == band]

def get_multi_band_photometry(observed_magnitudes):
    """Return a list of ObservedMagnitudes with duplicates averaged"""
    multi_band_photometry = []
    for obs in yield_observed_magnitudes_at_each_observed_band(observed_magnitudes):
        if len(obs) == 1:
            multi_band_photometry.append(obs[0])
        else:
            multi_band_photometry.append(combine_observed_magnitudes(obs))
    return multi_band_photometry

def get_observed_band_names(mutli_band_photometry_set):
    band_names = []
    for multi_band_photometry in mutli_band_photometry_set:
        band_names += [obs.band for obs in multi_band_photometry]
    return sorted(list(set(band_names)))

def get_observed_times(multi_band_photometry_set):
    times = []
    for multi_band_photometry in multi_band_photometry_set:
        times += get_times(multi_band_photometry)
    return sorted(list(set(times)))

def get_unobserved_times(lightcurve, observed_times):
    times = get_times(lightcurve)
    unobserved_times = []
    for observed_time in observed_times:
        if observed_time not in times:
            unobserved_times.append(observed_time)
    return unobserved_times

def get_times(observed_magnitudes):
    return [obs.time for obs in observed_magnitudes]

def get_lightcurve(multi_band_photometry_set, band):
    """Get list of all observed magnitudes in a given band """
    lightcurve = []
    for multi_band_photometry in multi_band_photometry_set:
        lightcurve += [obs for obs in multi_band_photometry if obs.band == band]
    return lightcurve
