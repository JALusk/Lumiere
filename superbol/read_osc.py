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

    magnitude = float(osc_photometry_dict['magnitude'])
    band = get_band(osc_photometry_dict['band'])
    
    if 'e_magnitude' in osc_photometry_dict.keys():
        uncertainty = float(osc_photometry_dict['e_magnitude'])
    else:
        uncertainty = 0.0

    return mag2flux.ObservedMagnitude(magnitude, uncertainty, band)

def get_band(band_name):
    """Make a Band object using the band_name given"""
    band_dict = retrieve_band_dict(band_name)
    effective_wavelength = band_dict['effective_wavelength']
    flux_conversion_factor = band_dict['flux_conversion_factor']
    return mag2flux.Band(band_name, 
                         effective_wavelength, 
                         flux_conversion_factor)

def retrieve_band_dict(band_name, path=None):
    """Load the Band attributes from a JSON file"""
    try:
        with open(path, 'r') as data_file:
            file_content = json.load(data_file)
        
        return file_content[band_name]
    except KeyError:
        raise NoBandFound

def retrieve_osc_photometry(sn_name, path=None):
    """Load the SN photometry from an OSC-formatted JSON file"""
    with open(path, 'r') as osc_data_file:
        file_content = json.load(osc_data_file)
    return file_content[sn_name]["photometry"]
