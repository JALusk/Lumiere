import math
import numpy as np

from itertools import groupby
from superbol import mag2flux
from operator import attrgetter

from scipy.interpolate import interp1d


class MissingFluxOutOfBounds(Exception):
    pass


def get_SED(fluxes):
    """Return a list of MonochromaticFluxes with duplicates averaged"""
    return sorted([combine_fluxes(flux) for flux in yield_fluxes_at_each_observed_wavelength(fluxes)], key=attrgetter("wavelength"))


def yield_fluxes_at_each_observed_wavelength(fluxes):
    """Yield lists of MonochromaticFluxes with the same wavelength"""
    for wavelength in set(f.wavelength for f in fluxes):
        yield [f for f in fluxes if f.wavelength == wavelength]


def combine_fluxes(fluxes):
    """Combine a list of MonochromaticFluxes using a weighted average"""
    values = get_flux_values(fluxes)
    uncertainties = get_flux_uncertainties(fluxes)
    combined_flux = weighted_average(values, uncertainties)
    combined_uncertainty = weighted_average_uncertainty(uncertainties)
    wavelength = fluxes[0].wavelength
    time = fluxes[0].time
    res = mag2flux.MonochromaticFlux(combined_flux,
                                     combined_uncertainty,
                                     wavelength,
                                     time)
    return res


def get_flux_values(fluxes):
    return [f.flux for f in fluxes]


# TODO Create test for passing flux with uncertainty of 0
def get_flux_uncertainties(fluxes):
    return [f.flux_uncertainty for f in fluxes]


def get_times(fluxes):
    return [f.time for f in fluxes]


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
    return [1/s**2 for s in uncertainties]


def group_fluxes(fluxes, keyfunc=math.floor):
    """Group fluxes by applying keyfunc() to each flux.time in the fluxes.

    Fluxes with the same return value for keyfunc(flux.time) will be 
    part of the same group."""
    grouped_fluxes = [list(it) for k, it in groupby(
        fluxes, lambda x: keyfunc(x.time))]

    for flux_group in grouped_fluxes:
        for flux in flux_group:
            flux.time = keyfunc(flux.time)

    return grouped_fluxes


def interpolate_missing_fluxes(SEDs):
    observed_wavelengths = get_observed_wavelengths(SEDs)
    observed_times = get_observed_times(SEDs)
    for wl in observed_wavelengths:
        monochromatic_lightcurve = get_monochromatic_lightcurve(SEDs, wl)
        interpolated_fluxes = get_interpolated_fluxes(
            monochromatic_lightcurve, observed_times)
        append_interpolated_fluxes_to_SEDs(interpolated_fluxes, SEDs)


def append_interpolated_fluxes_to_SEDs(interpolated_fluxes, SEDs):
    for SED in SEDs:
        for interpolated_flux in interpolated_fluxes:
            if interpolated_flux.time == SED[0].time:
                SED.append(interpolated_flux)


def get_interpolated_fluxes(lightcurve, observed_times):
    unobserved_times = get_unobserved_times(lightcurve, observed_times)
    interpolated_fluxes = []
    for unobserved_time in unobserved_times:
        try:
            previous_flux = get_previous_flux(lightcurve, unobserved_time)
            next_flux = get_next_flux(lightcurve, unobserved_time)
            if next_flux.time - previous_flux.time <= 2:
                f = interp1d([previous_flux.time, next_flux.time],
                             [previous_flux.flux, next_flux.flux])
                interpolated_flux_value = f(unobserved_time)
                interpolated_flux_uncertainty = get_interpolated_flux_uncertainty(
                    previous_flux,
                    next_flux,
                    unobserved_time)
                interpolated_flux = mag2flux.MonochromaticFlux(
                    interpolated_flux_value, interpolated_flux_uncertainty, previous_flux.wavelength, unobserved_time)
                interpolated_fluxes.append(interpolated_flux)
        except MissingFluxOutOfBounds:
            pass
    return interpolated_fluxes


def get_interpolated_flux_uncertainty(previous_flux, next_flux, unobserved_time):
    weight1 = (next_flux.time - unobserved_time) / \
        (next_flux.time - previous_flux.time)
    weight2 = (unobserved_time - previous_flux.time) / \
        (next_flux.time - previous_flux.time)
    return math.sqrt(weight1**2 * previous_flux.flux_uncertainty**2 + weight2**2 + next_flux.flux_uncertainty**2)


# TODO Code review
def get_previous_flux(monochromatic_lightcurve, unobserved_time):
    earlier_fluxes = []
    for flux in monochromatic_lightcurve:
        delta = unobserved_time - flux.time
        if delta > 0:
            earlier_fluxes.append(flux)
    if earlier_fluxes == []:
        raise MissingFluxOutOfBounds
    return max(earlier_fluxes, key=lambda flux: flux.time)


def get_next_flux(monochromatic_lightcurve, unobserved_time):
    later_fluxes = [
        flux for flux in monochromatic_lightcurve if unobserved_time - flux.time < 0]
    if later_fluxes == []:
        raise MissingFluxOutOfBounds
    return min(later_fluxes, key=lambda flux: flux.time)


# TODO Code review
def get_gap_size(times, unobserved_time):
    positive_deltas = [
        time - unobserved_time for time in times if time - unobserved_time > 0]
    negative_deltas = [
        time - unobserved_time for time in times if time - unobserved_time <= 0]
    return min(positive_deltas) - max(negative_deltas)


def get_unobserved_times(lightcurve, observed_times):
    times = get_times(lightcurve)
    return [
        observed_time for observed_time in observed_times if observed_time not in times]


def get_observed_wavelengths(SEDs):
    wavelengths = []
    for SED in SEDs:
        wavelengths += [f.wavelength for f in SED]
    return sorted(list(set(wavelengths)))


def get_observed_times(SEDs):
    times = []
    for SED in SEDs:
        times += get_times(SED)
    return sorted(list(set(times)))


def get_monochromatic_lightcurve(SEDs, wavelength):
    """Get list of all fluxes at a given wavelength across all SEDs"""
    monochromatic_lightcurve = []
    for SED in SEDs:
        monochromatic_lightcurve += [
            f for f in SED if f.wavelength == wavelength]
    return monochromatic_lightcurve
