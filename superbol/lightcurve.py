import math

from superbol import sed
from superbol import photometry

# TODO Doc these
# TODO No tests written for fns in this file
def calculate_lightcurve(fluxes, distance, flux_calculator):
    grouped_fluxes = sed.group_fluxes(fluxes, math.floor)
    SEDs = [sed.get_SED(flux_group)
            for flux_group in grouped_fluxes if len(flux_group) > 2]

    sed.interpolate_missing_fluxes(SEDs)

    return [flux_calculator(SED).to_lbol(distance) for SED in SEDs]


def calculate_bc_lightcurve(magnitudes, distance, flux_calculator):
    grouped_magnitudes = photometry.group_magnitudes(magnitudes, math.floor)
    multi_band_photometry_set = [
        photometry.get_multi_band_photometry(magnitude_group) for magnitude_group in grouped_magnitudes if len(magnitude_group) > 2]

    photometry.interpolate_missing_magnitudes(multi_band_photometry_set)
    return [flux_calculator(multi_band_photometry).to_lbol(distance) for multi_band_photometry in multi_band_photometry_set]
