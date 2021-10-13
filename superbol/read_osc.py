import json

from superbol import mag2flux

class NoMagnitude(Exception):
    pass

class NoBandNameGiven(Exception):
    pass

class NoTimeGiven(Exception):
    pass

class NoBandFound(Exception):
    pass

class NoUncertainty(Exception):
    pass

def get_observed_magnitude(osc_photometry_dict):
    """Turn photometry from the OSC into an ObservedMagnitude"""
    try:
        magnitude = float(osc_photometry_dict['magnitude'])
        uncertainty = float(osc_photometry_dict['e_magnitude'])
        band = get_band(osc_photometry_dict['band'])
        time = float(osc_photometry_dict['time'])
        source = osc_photometry_dict['source']
    except KeyError:
        if 'magnitude' not in osc_photometry_dict.keys():
            raise NoMagnitude
        if 'e_magnitude' not in osc_photometry_dict.keys():
            raise NoUncertainty
        elif 'band' not in osc_photometry_dict.keys():
            raise NoBandNameGiven
        elif 'time' not in osc_photometry_dict.keys():
            raise NoTimeGiven

    observed_magnitude = mag2flux.ObservedMagnitude(magnitude, uncertainty, band, time)
    observed_magnitude.source = source

    return observed_magnitude

def get_band(band_name):
    """Make a Band object using the band_name given"""
    band_dict = retrieve_band_dict(band_name)
    band_alt_name = band_dict['alt_name']
    effective_wavelength = band_dict['effective_wavelength']
    flux_conversion_factor = band_dict['flux_conversion_factor']
    return mag2flux.Band(band_name,
                         band_alt_name,
                         effective_wavelength, 
                         flux_conversion_factor)

def retrieve_band_dict(band_name, path='./data/bands.json'):
    """Load the Band attributes from a JSON file"""
    try:
        with open(path, 'r') as data_file:
            file_content = json.load(data_file)
        
        return file_content[band_name]
    except KeyError:
        raise NoBandFound("The band {0} was not found in the SuperBoL band catalog.".format(band_name))

# TODO Add default value?
def retrieve_osc_photometry(sn_name, path=None):
    """Load the SN photometry from an OSC-formatted JSON file"""
    with open(path, 'r') as osc_data_file:
        file_content = json.load(osc_data_file)
    return file_content[sn_name]["photometry"]
