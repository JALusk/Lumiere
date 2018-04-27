import math

from superbol import sed
from superbol import photometry

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

def calculate_bc_lightcurve(magnitudes, distance, luminosity_calculator):
    grouped_magnitudes = photometry.group_magnitudes(magnitudes, math.floor)
    lightcurve = []
    multi_band_photometry_set = []
    for magnitude_group in grouped_magnitudes:
        if len(magnitude_group) > 2:
            multi_band_photometry_set.append(photometry.get_multi_band_photometry(magnitude_group))

    photometry.interpolate_missing_magnitudes(multi_band_photometry_set)

    for multi_band_photometry in multi_band_photometry_set:
        luminosity = luminosity_calculator(multi_band_photometry, distance)
        lightcurve.append(luminosity)
    return lightcurve
