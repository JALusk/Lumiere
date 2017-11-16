import math
import numpy as np

from superbol import mag2flux

class InsufficientFluxes(Exception):
    pass

class Distance(object):

    def __init__(self, value, uncertainty):
        self.value = value
        self.uncertainty = uncertainty

class QuasiBolometricFlux(object):

    def __init__(self, value, uncertainty, time):
        self.value = value
        self.uncertainty = uncertainty
        self.time = time

class QuasiBolometricLuminosity(object):

    def __init__(self, value, uncertainty, time):
        self.value = value
        self.uncertainty = uncertainty
        self.time = time

def combine_fluxes(fluxes):
    """Combine a list of MonochromaticFluxes into one average MonochromaticFlux"""
    combined_flux = np.mean([f.flux for f in fluxes])
    combined_uncertainty = np.sqrt(sum(f.flux_uncertainty**2 for f in fluxes)) / len(fluxes)
    wavelength = fluxes[0].wavelength
    time = fluxes[0].time
    return mag2flux.MonochromaticFlux(combined_flux,
                                      combined_uncertainty,
                                      wavelength,
                                      time)

def yield_fluxes_at_each_observed_wavelength(fluxes):
    """Yield lists of MonochromaticFluxes with the same wavelength"""
    for wavelength in set(f.wavelength for f in fluxes):
        yield [f for f in fluxes if f.wavelength == wavelength]

def get_integrable_fluxes(fluxes):
    """Return a list of MonochromaticFluxes with duplicates averaged"""
    integrable_fluxes = []
    for f in yield_fluxes_at_each_observed_wavelength(fluxes):
        if len(f) == 1:
            integrable_fluxes.append(f[0])
        else:
            integrable_fluxes.append(combine_fluxes(f))
    return integrable_fluxes

def get_quasi_bolometric_flux(integral_calculator, 
                              uncertainty_calculator, 
                              fluxes):
    """Calculate Fqbol using the supplied integration technique"""
    integrable_fluxes = get_integrable_fluxes(fluxes)
    if len(integrable_fluxes) < 2:
        raise InsufficientFluxes(
            "Cannot calculate quasi-bolometric flux with fewer " +
            "than two fluxes, {0} received".format(len(fluxes)))
    
    return QuasiBolometricFlux(
            value = integral_calculator.calculate(integrable_fluxes),
            uncertainty = uncertainty_calculator(integrable_fluxes),
            time = fluxes[0].time)

class TrapezoidalIntegralCalculator(object):
    """Integrate between fluxes using the trapezoidal method"""
    def calculate(self, fluxes):
        """Calculate the integral using numpy.trapz"""
        self._sort_fluxes_by_wavelength(fluxes)
        flux_list = self._get_flux_list(fluxes)
        wavelength_list = self._get_wavelength_list(fluxes)
        return np.trapz(flux_list, wavelength_list)

    def _sort_fluxes_by_wavelength(self, fluxes):
        """Sort the fluxes in-place by wavelength"""
        fluxes.sort(key = lambda x: x.wavelength)

    def _get_flux_list(self, fluxes):
        """Return a list of flux values"""
        return [f.flux for f in fluxes]

    def _get_wavelength_list(self, fluxes):
        """Return a list of flux wavelengths"""
        return [f.wavelength for f in fluxes]

def uncertainty_calculator_trapezoidal(fluxes):
    """Calculate uncertainty in trapezoidal integral of fluxes"""
    radicand = 0

    for i, flux in enumerate(fluxes):
        if i == 0:
            radicand += (0.5 * (fluxes[i+1].wavelength - flux.wavelength) 
                         * flux.flux_uncertainty)**2
        elif i == len(fluxes) - 1:
            radicand += (0.5 * (flux.wavelength - fluxes[i-1].wavelength)
                         * flux.flux_uncertainty)**2
        else:
            radicand += (0.5 * (fluxes[i+1].wavelength - fluxes[i-1].wavelength)
                         * flux.flux_uncertainty)**2

    return math.sqrt(radicand)

def convert_flux_to_luminosity(fqbol, distance):
    """Convert quasi-bolometric flux to quasi-bolometric luminosity"""
    lbqol_value = fqbol.value * 4.0 * math.pi * distance.value**2
    lqbol_uncertainty = math.sqrt((4.0 * math.pi * distance.value**2 * fqbol.uncertainty)**2 + (fqbol.value * 8.0 * math.pi * distance.value * distance.uncertainty)**2)

    return QuasiBolometricLuminosity(lbqol_value, lqbol_uncertainty, fqbol.time)

def calculate_qbol_luminosity(flux_group, distance):
    """Turn a group of fluxes into a quasi-bolometric luminosity"""
    integral_calculator = TrapezoidalIntegralCalculator()
    uncertainty_calculator = uncertainty_calculator_trapezoidal
    
    fqbol = get_quasi_bolometric_flux(integral_calculator,
                                      uncertainty_calculator,
                                      flux_group)
    lqbol = convert_flux_to_luminosity(fqbol, distance)
    
    return lqbol
