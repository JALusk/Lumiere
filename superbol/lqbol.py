import math
import numpy as np

class InsufficientFluxes(Exception):
    pass

class Distance(object):

    def __init__(self, value, uncertainty):
        self.value = value
        self.uncertainty = uncertainty

class QuasiBolometricFlux(object):

    def __init__(self, value, uncertainty):
        self.value = value
        self.uncertainty = uncertainty

class QuasiBolometricLuminosity(object):

    def __init__(self, value, uncertainty):
        self.value = value
        self.uncertainty = uncertainty

def get_quasi_bolometric_flux(integral_calculator, 
                              uncertainty_calculator, 
                              fluxes):
    """Calculate Fqbol using the supplied integration technique"""
    if len(fluxes) < 2:
        raise InsufficientFluxes(
            "Cannot calculate quasi-bolometric flux with fewer " +
            "than two fluxes, {0} received".format(len(fluxes)))

    return QuasiBolometricFlux(
            value = integral_calculator.calculate(fluxes),
            uncertainty = uncertainty_calculator(fluxes))

class TrapezoidalIntegralCalculator(object):
    """Integrate between fluxes using the trapezoidal method"""
    def _sort_fluxes_by_wavelength(self, fluxes):
        """Sort the fluxes in-place by wavelength"""
        fluxes.sort(key = lambda x: x.wavelength)

    def _get_flux_list(self, fluxes):
        """Return a list of flux values"""
        return [f.flux for f in fluxes]

    def _get_wavelength_list(self, fluxes):
        """Return a list of flux wavelengths"""
        return [f.wavelength for f in fluxes]

    def calculate(self, fluxes):
        """Calculate the integral using numpy.trapz"""
        self._sort_fluxes_by_wavelength(fluxes)
        flux_list = self._get_flux_list(fluxes)
        wavelength_list = self._get_wavelength_list(fluxes)
        return np.trapz(flux_list, wavelength_list)

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

    return QuasiBolometricLuminosity(lbqol_value, lqbol_uncertainty)

def calculate_qbol_luminosity(flux_group, distance):
    """Turn a group of fluxes into a quasi-bolometric luminosity"""
    integral_calculator = TrapezoidalIntegralCalculator()
    uncertainty_calculator = uncertainty_calculator_trapezoidal
    
    fqbol = get_quasi_bolometric_flux(integral_calculator,
                                      uncertainty_calculator,
                                      flux_group)
    lqbol = convert_flux_to_luminosity(fqbol, distance)
    
    return lqbol
