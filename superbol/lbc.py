import math
import json
import numpy as np
import itertools

from astropy import units as u

class BCColorRelation(object):

    def get_coefficients(self, band1, band2):
        """Get coefficients for color combination band1 - band2"""
        color = band1.name + "-" + band2.name
        if color == "B-V":
            return self.coefficients_BV()
        elif color == "V-I":
            return self.coefficients_VI()

    def coefficients_BV(self):
        return [-0.823, 5.027, -13.409, 20.133, -18.096, 9.084, -1.950]

    def coefficients_VI(self):
        return [-1.355, 6.262, -2.676, -22.973, 35.542, -15.340]

def calculate_bc_luminosity(obs_group, distance):
    """Turn a group of observations into an average BC luminosity"""
    obs_group.sort(key=lambda x: x.band.effective_wavelength)
    if 'V' in [x.band.name for x in obs_group]:
        lbc = []
        for obs_pair in itertools.combinations(obs_group, 2):
            try:
                bc = compute_bolometric_correction('BH09', obs_pair[0], obs_pair[1])
                mbol = apply_bolometric_correction(bc, next(obs for obs in obs_group if obs.band.name == 'V'))
                Mbol = convert_apparent_to_absolute_magnitude(mbol, distance)
                lbc.append(convert_Mbol_to_Lbol(Mbol))
            except InvalidColor:
                pass
            except InvalidBCMethod:
                pass
            except InvalidFilterCombination:
                pass

def compute_bolometric_correction(method_name, obs1, obs2):
    """Calculate the bolometric correction for obs1-obs2 using specified method"""
    bc_method_data = retrieve_bc_method_data(method_name)
    color = obs1.magnitude - obs2.magnitude
    coefficients = get_bc_method_coefficients(bc_method_data, obs1.band, obs2.band)
    range_min, range_max = get_bc_method_range(bc_method_data, obs1.band, obs2.band)
    rms = get_bc_method_rms(bc_method_data, obs1.band, obs2.band)
    if range_min <= color <= range_max:
        bc = compute_polynomial(color, coefficients)
        uncertainty = rms
        return BolometricCorrection(bc, rms)
    else:
        raise InvalidColor("Cannot calculate bolometric correction, given " +
                           "{0}-{1} color outside the allowable range for " +
                           "the {2} method.".format(obs1.band.name, 
                                                    obs2.band.name,
                                                    method_name))

def retrieve_bc_method_data(method_name, path="/home/jlusk/src/superbol/data/bc_color_fits.json"):
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


def apply_bolometric_correction(bc, observed_magnitude):
    """Apply the bolometric correction to the observed magnitude"""
    mbol_value = bc.value + observed_magnitude.magnitude
    mbol_uncertainty = math.sqrt(bc.uncertainty**2 + observed_magnitude.uncertainty**2)
    mbol_time = observed_magnitude.time
    return BolometricMagnitude(mbol_value, mbol_uncertainty, mbol_time)


    return BCLuminosity(np.mean([x.value for x in lbc]), np.mean([x.uncertainty for x in lbc]), np.mean([x.time for x in lbc]))

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

class BCLuminosity(object):

    def __init__(self, value, uncertainty, time):
        self.value = value
        self.uncertainty = uncertainty
        self.time = time

class BolometricCorrection(object):

    def __init__(self, value, uncertainty):
        self.value = value
        self.uncertainty = uncertainty

class AbsoluteMagnitude(object):

    def __init__(self, value, uncertainty, time):
        self.value = value
        self.uncertainty = uncertainty
        self.time = time


def compute_polynomial(color, coefficients):
    """Compute the bc-color polynomial"""
    result = 0.0
    for n, coefficient in enumerate(coefficients):
        result += coefficient * color**n
    return result


def convert_Mbol_to_Lbol(bolometric_magnitude):
    L0 = 3.0128E35
    L = L0 * 10**(-0.4 * bolometric_magnitude.value)
    L_uncertainty = -0.4 * math.log(10) * L * bolometric_magnitude.uncertainty
    return BCLuminosity(L, L_uncertainty, bolometric_magnitude.time)

def convert_apparent_to_absolute_magnitude(apparent_magnitude, distance):
    distance_pc = distance.value / 3.086E+18
    absolute_magnitude = apparent_magnitude.value - 5.0 * math.log(distance_pc / 10.0, 10)
    uncertainty = math.sqrt(apparent_magnitude.uncertainty**2 + (5 / (distance_pc * math.log(10)) * distance.uncertainty / 3.086E+18)**2)
    return AbsoluteMagnitude(absolute_magnitude, uncertainty, apparent_magnitude.time)
