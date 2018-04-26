import math
import numpy as np

from itertools import groupby
from superbol import mag2flux

from scipy.interpolate import interp1d

class MissingMagnitudeOutOfBounds(Exception):
    pass

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

def interpolate_missing_magnitudes(multi_band_photometry_set):
    observed_bands = get_observed_band_names(multi_band_photometry_set)
    observed_times = get_observed_times(multi_band_photometry_set)
    for band in observed_bands:
        lightcurve = get_lightcurve(multi_band_photometry_set, band)
        interpolated_magnitudes = get_interpolated_magnitudes(lightcurve, observed_times)
        append_interpolated_magnitudes_to_multi_band_photometry_set(interpolated_magnitudes,
                multi_band_photometry_set)

def get_interpolated_magnitudes(lightcurve, observed_times):
    unobserved_times = get_unobserved_times(lightcurve, observed_times)
    interpolated_magnitudes = []
    for unobserved_time in unobserved_times:
        try:
            previous_observed_magnitude = get_previous_observed_magnitude(lightcurve, unobserved_time)
            next_observed_magnitude = get_next_observed_magnitude(lightcurve, unobserved_time)
            if next_observed_magnitude.time - previous_observed_magnitude.time <= 2:
                f = interp1d([previous_observed_magnitude.time, next_observed_magnitude.time], 
                          [previous_observed_magnitude.magnitude, next_observed_magnitude.magnitude])
                interpolated_magnitude_value = f(unobserved_time)
                interpolated_magnitude_uncertainty = get_interpolated_magnitude_uncertainty(
                                                 previous_observed_magnitude,
                                                 next_observed_magnitude,
                                                 unobserved_time)
                interpolated_magnitude = mag2flux.ObservedMagnitude(interpolated_magnitude_value, interpolated_magnitude_uncertainty, previous_observed_magnitude.band, unobserved_time)
                interpolated_magnitudes.append(interpolated_magnitude)
        except MissingMagnitudeOutOfBounds:
            pass
    return interpolated_magnitudes

def get_interpolated_magnitude_uncertainty(previous_observed_magnitude, next_observed_magnitude, unobserved_time):
    weight1 = (next_observed_magnitude.time - unobserved_time)/(next_observed_magnitude.time - previous_observed_magnitude.time)
    weight2 = (unobserved_time - previous_observed_magnitude.time)/(next_observed_magnitude.time - previous_observed_magnitude.time)
    return math.sqrt(weight1**2 * previous_observed_magnitude.uncertainty**2 + weight2**2 + next_observed_magnitude.uncertainty**2)

def append_interpolated_magnitudes_to_multi_band_photometry_set(interpolated_magnitudes, multi_band_photometry_set):
    for multi_band_photometry in multi_band_photometry_set:
        for interpolated_magnitude in interpolated_magnitudes:
            if interpolated_magnitude.time == multi_band_photometry[0].time:
                multi_band_photometry.append(interpolated_magnitude)

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

def get_previous_observed_magnitude(lightcurve, unobserved_time):
    """Return the most recently observed magnitude before the unobserved time"""
    earlier_observed_magnitudes = []
    for observed_magnitude in lightcurve:
        delta = unobserved_time - observed_magnitude.time
        if delta > 0:
            earlier_observed_magnitudes.append(observed_magnitude)
    if earlier_observed_magnitudes == []:
        raise MissingMagnitudeOutOfBounds
    return max(earlier_observed_magnitudes, key=lambda observed_magnitude: observed_magnitude.time)

def get_next_observed_magnitude(lightcurve, unobserved_time):
    """Return the most recently observed magnitude after the unobserved time"""
    later_observed_magnitudes = []
    for observed_magnitude in lightcurve:
        delta = unobserved_time - observed_magnitude.time
        if delta < 0:
            later_observed_magnitudes.append(observed_magnitude)
    if later_observed_magnitudes == []:
        raise MissingMagnitudeOutOfBounds
    return min(later_observed_magnitudes, key=lambda observed_magnitude: observed_magnitude.time)

