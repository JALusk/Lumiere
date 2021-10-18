import math
import itertools
import json
import numpy as np
import os

from astropy import units as u
from superbol.lum import BolometricFlux

dirname = os.path.dirname(__file__)
bc_color_fits = os.path.join(dirname, '../data/bc_color_fits.json')

def calculate_bc_flux_h01(obs_group):
    """Turn a group of observations into an average BC flux
    
    Args: 
        obs_group: observed photometry

    Returns: 
        BC Bolemetric flux
    """  
    obs_group.sort(key=lambda x: x.band.effective_wavelength)
    if 'V' in [x.band.name for x in obs_group]:
        fbc = []
        for obs_pair in itertools.combinations(obs_group, 2):
            try:
                bc = compute_bolometric_correction('H01', obs_pair[0], obs_pair[1])
                mbol = apply_bolometric_correction(bc, next(obs for obs in obs_group if obs.band.name == 'V'))
                zeropoint = -10.88802466
                Fbol = convert_mbol_to_Fbol(mbol, zeropoint)
                fbc.append(Fbol)
            except InvalidColor:
                pass
            except InvalidBCMethod:
                pass
            except InvalidFilterCombination:
                pass
    return BCBolometricFlux(np.mean([x.value for x in fbc]), np.mean([x.uncertainty for x in fbc]), np.mean([x.time for x in fbc]))

def calculate_bc_flux_bh09(obs_group):
    """Turn a group of observations into an average BC flux"""
    obs_group.sort(key=lambda x: x.band.effective_wavelength)
    if 'V' in [x.band.name for x in obs_group]:
        fbc = []
        for obs_pair in itertools.combinations(obs_group, 2):
            try:
                bc = compute_bolometric_correction('BH09', obs_pair[0], obs_pair[1])
                mbol = apply_bolometric_correction(bc, next(obs for obs in obs_group if obs.band.name == 'V'))
                zeropoint = -11.64
                Fbol = convert_mbol_to_Fbol(mbol, zeropoint)
                fbc.append(Fbol)
            except InvalidColor:
                pass
            except InvalidBCMethod:
                pass
            except InvalidFilterCombination:
                pass
    return BCBolometricFlux(np.mean([x.value for x in fbc]), np.mean([x.uncertainty for x in fbc]), np.mean([x.time for x in fbc]))


def compute_bolometric_correction(method_name, obs1, obs2):
    """Calculate the bolometric correction for obs1-obs2 using specified method"""
    bc_method_data = retrieve_bc_method_data(method_name)
    color = obs1.magnitude - obs2.magnitude
    color_err = np.sqrt(obs1.uncertainty**2 + obs2.uncertainty**2)
    coefficients = get_bc_method_coefficients(bc_method_data, obs1.band, obs2.band)
    range_min, range_max = get_bc_method_range(bc_method_data, obs1.band, obs2.band)
    rms = get_bc_method_rms(bc_method_data, obs1.band, obs2.band)
    if range_min <= color <= range_max:
        bc = compute_polynomial(color, coefficients)
        poly_deriv = compute_polynomial_derivative(color, coefficients)
        uncertainty = np.sqrt(rms**2 + (poly_deriv * color_err)**2)
        return BolometricCorrection(bc, uncertainty)
    else:
        raise InvalidColor("Cannot calculate bolometric correction, given " +
                           "{0}-{1} color outside the allowable range for " +
                           "the {2} method.".format(obs1.band.name, 
                                                    obs2.band.name,
                                                    method_name))

def retrieve_bc_method_data(method_name, path=bc_color_fits):
    with open(path, 'r') as bc_color_data_file:
        file_content = json.load(bc_color_data_file)
        try:
            return file_content[method_name]
        except KeyError:
            raise InvalidBCMethod("Bolometric correction method {0} not found".format(method_name))

def get_bc_method_coefficients(method_data, band1, band2):
    try:
        return method_data[band1.name + '-' + band2.name]['coefficients']
    except KeyError:
        raise InvalidFilterCombination("{0} - {1} not an allowable color for this method".format(band1.name, band2.name))

def get_bc_method_range(method_data, band1, band2):
    range_min = method_data[band1.name + '-' + band2.name]['range_min']
    range_max = method_data[band1.name + '-' + band2.name]['range_max']
    return range_min, range_max

def get_bc_method_rms(method_data, band1, band2):
    return method_data[band1.name + '-' + band2.name]['rms']

def get_bc_method_zeropoint(method_data):
    return method_data['properties']['ZP']

def apply_bolometric_correction(bc, observed_magnitude):
    """Apply the bolometric correction to the observed magnitude"""
    mbol_value = bc.value + observed_magnitude.magnitude
    mbol_uncertainty = math.sqrt(bc.uncertainty**2 + observed_magnitude.uncertainty**2)
    mbol_time = observed_magnitude.time
    return BolometricMagnitude(mbol_value, mbol_uncertainty, mbol_time)


class InvalidColor(Exception):
    pass

class InvalidBCMethod(Exception):
    pass

class InvalidFilterCombination(Exception):
    pass

class BolometricMagnitude(object):

    def __init__(self, value, uncertainty, time):
        self.value = value
        self.uncertainty = uncertainty
        self.time = time

class BCBolometricFlux(BolometricFlux):

    def __init__(self, value, uncertainty, time):
        super().__init__(value, uncertainty)
        self.time = time
   
    def to_lbol(self, distance):
        lbol = super().to_lbol(distance)
        lbol.time = self.time
        return lbol


class BolometricCorrection(object):

    def __init__(self, value, uncertainty):
        self.value = value
        self.uncertainty = uncertainty

def compute_polynomial(color, coefficients):
    """Compute the bc-color polynomial"""
    result = 0.0
    for n, coefficient in enumerate(coefficients):
        result += coefficient * color**n
    return result

def compute_polynomial_derivative(color, coefficients):
    """Compute the derivative of the bc-color polynomial"""
    result = 0.0
    if color == 0.0:
        result = coefficients[1]
    else:    
        for n, coefficient in enumerate(coefficients):
            result += n * coefficient * color**(n-1)
    return result

def convert_mbol_to_Fbol(mbol, zeropoint):
    Fbol = 10**((-mbol.value + zeropoint)/2.5)
    Fbol_uncertainty = np.abs(math.log(10)/2.5 * Fbol * mbol.uncertainty)
    return BCBolometricFlux(Fbol, Fbol_uncertainty, mbol.time)
