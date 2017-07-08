import json

from superbol import mag2flux

class NoMagnitude(Exception):
    pass

class NoBandNameGiven(Exception):
    pass

class NoBandFound(Exception):
    pass

def get_observed_magnitude(osc_photometry_dict):
    """Turn photometry from the OSC into an ObservedMagnitude"""
    if 'magnitude' not in osc_photometry_dict.keys():
        raise NoMagnitude
    elif 'band' not in osc_photometry_dict.keys():
        raise NoBandNameGiven

    magnitude = osc_photometry_dict['magnitude']
    uncertainty = osc_photometry_dict['e_magnitude']
    band = get_band(osc_photometry_dict['band'])
    return mag2flux.ObservedMagnitude(magnitude, uncertainty, band)

def get_band(band_name):
    band_dict = retrieve_band_dict(band_name)
    effective_wavelength = band_dict['effective_wavelength']
    flux_conversion_factor = band_dict['flux_conversion_factor']
    return mag2flux.Band(band_name, 
                         effective_wavelength, 
                         flux_conversion_factor)

def retrieve_band_dict(band_name, path=None):
    try:
        with open(path, 'r') as data_file:
            file_content = json.load(data_file)
        
        return file_content[band_name]
    except KeyError:
        raise NoBandFound
