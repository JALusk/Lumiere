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
    return [list(it) for k, it in groupby(
        magnitudes, lambda x: keyfunc(x.time))]


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
    return 1.0/math.sqrt(sum(weights))


def get_weights(uncertainties):
    """Calculate weights from uncertainties"""
    return [1/s**2 for s in uncertainties]


def yield_observed_magnitudes_at_each_observed_band(observed_magnitudes):
    """Yield lists of ObservedMagnitudes with the same band"""
    list_of_band_names = sorted(
        list(set(obs.band.name for obs in observed_magnitudes)))
    for band_name in list_of_band_names:
        yield [obs for obs in observed_magnitudes if obs.band.name == band_name]


def get_multi_band_photometry(observed_magnitudes):
    """Return a list of ObservedMagnitudes with duplicates averaged"""
    return [
        obs[0] if len(obs) == 1 else combine_observed_magnitudes(obs)
        for obs in yield_observed_magnitudes_at_each_observed_band(observed_magnitudes)
    ]


def interpolate_missing_magnitudes(multi_band_photometry_set):
    observed_bands = get_observed_band_names(multi_band_photometry_set)
    observed_times = get_observed_times(multi_band_photometry_set)
    for band in observed_bands:
        lightcurve = get_lightcurve(multi_band_photometry_set, band)
        interpolated_magnitudes = get_interpolated_magnitudes(
            lightcurve, observed_times)
        append_interpolated_magnitudes_to_multi_band_photometry_set(interpolated_magnitudes,
                                                                    multi_band_photometry_set)


def get_interpolated_magnitudes(lightcurve, observed_times):
    unobserved_times = get_unobserved_times(lightcurve, observed_times)
    interpolated_magnitudes = []
    for unobserved_time in unobserved_times:
        try:
            previous_observed_magnitude = get_previous_observed_magnitude(
                lightcurve, unobserved_time)
            next_observed_magnitude = get_next_observed_magnitude(
                lightcurve, unobserved_time)
            if next_observed_magnitude.time - previous_observed_magnitude.time <= 2:
                f = interp1d([previous_observed_magnitude.time, next_observed_magnitude.time],
                             [previous_observed_magnitude.magnitude, next_observed_magnitude.magnitude])
                interpolated_magnitude_value = f(unobserved_time)
                interpolated_magnitude_uncertainty = get_interpolated_magnitude_uncertainty(
                    previous_observed_magnitude,
                    next_observed_magnitude,
                    unobserved_time)
                interpolated_magnitude = mag2flux.ObservedMagnitude(
                    interpolated_magnitude_value, interpolated_magnitude_uncertainty, previous_observed_magnitude.band, unobserved_time)
                interpolated_magnitudes.append(interpolated_magnitude)
        except MissingMagnitudeOutOfBounds:
            # TODO No test written
            pass
    return interpolated_magnitudes


def get_interpolated_magnitude_uncertainty(previous_observed_magnitude, next_observed_magnitude, unobserved_time):
    weight1 = (next_observed_magnitude.time - unobserved_time) / \
        (next_observed_magnitude.time - previous_observed_magnitude.time)
    weight2 = (unobserved_time - previous_observed_magnitude.time) / \
        (next_observed_magnitude.time - previous_observed_magnitude.time)
    return math.sqrt(weight1**2 * previous_observed_magnitude.uncertainty**2 + weight2**2 + next_observed_magnitude.uncertainty**2)


def append_interpolated_magnitudes_to_multi_band_photometry_set(interpolated_magnitudes, multi_band_photometry_set):
    for multi_band_photometry in multi_band_photometry_set:
        for interpolated_magnitude in interpolated_magnitudes:
            if interpolated_magnitude.time == multi_band_photometry[0].time:
                multi_band_photometry.append(interpolated_magnitude)


def get_observed_band_names(mutli_band_photometry_set):
    band_names = []
    for multi_band_photometry in mutli_band_photometry_set:
        band_names += [obs.band.name for obs in multi_band_photometry]
    return sorted(list(set(band_names)))


def get_observed_times(multi_band_photometry_set):
    times = []
    for multi_band_photometry in multi_band_photometry_set:
        times += get_times(multi_band_photometry)
    return sorted(list(set(times)))


def get_unobserved_times(lightcurve, observed_times):
    times = get_times(lightcurve)
    return [obs for obs in observed_times if obs not in times]


def get_times(observed_magnitudes):
    return [obs.time for obs in observed_magnitudes]


def get_lightcurve(multi_band_photometry_set, band_name):
    """Get list of all observed magnitudes in a given band """
    # TODO Find out why this filters band_name, maybe create a convenience function
    return [
        obs for multiband_photometry in multi_band_photometry_set for obs in multiband_photometry if obs.band.name == band_name
    ]


def get_previous_observed_magnitude(lightcurve, unobserved_time):
    """Return the most recently observed magnitude before the unobserved time"""
    earlier_observed_magnitudes = [
        obs for obs in lightcurve if unobserved_time - obs.time > 0]
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
