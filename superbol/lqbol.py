import math
import numpy as np

class InsufficientFluxes(Exception):
    pass

class QuasiBolometricFlux(object):
    def __init__(self, value, uncertainty):
        self.value = value
        self.uncertainty = uncertainty

class TrapezoidalIntegralCalculator(object):
    """Integrate between fluxes using the trapezoidal method"""
    def _get_flux_list(self, fluxes):
        """Return a list of flux values"""
        return [f.flux for f in fluxes]

    def _get_wavelength_list(self, fluxes):
        """Return a list of flux wavelengths"""
        return [f.wavelength for f in fluxes]

    def calculate(self, fluxes):
        """Calculate the integral using numpy.trapz"""
        flux_list = self._get_flux_list(fluxes)
        wavelength_list = self._get_wavelength_list(fluxes)
        return np.trapz(flux_list, wavelength_list)

def get_quasi_bolometric_flux(integral_calculator, 
                              uncertainty_calculator, 
                              fluxes):
    if len(fluxes) < 2:
        raise InsufficientFluxes(
            "Cannot calculate quasi-bolometric flux with fewer " +
            "than two fluxes, {0} received".format(len(fluxes)))

    return QuasiBolometricFlux(
            value = integral_calculator.calculate(fluxes),
            uncertainty = uncertainty_calculator(fluxes))

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
