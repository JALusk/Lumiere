import math
import numpy as np

from superbol import mag2flux
from superbol.lum import BolometricFlux

class InsufficientFluxes(Exception):
    pass

class QuasiBolometricFlux(BolometricFlux):

    def __init__(self, value, uncertainty, time):
        super().__init__(value, uncertainty)
        self.time = time

    def to_lbol(self, distance):
        lbol = super().to_lbol(distance)
        lbol.time = self.time
        return lbol

def get_quasi_bolometric_flux(integral_calculator, 
                              uncertainty_calculator, 
                              SED):
    """Calculate Fqbol using the supplied integration technique"""
    if len(SED) < 2:
        raise InsufficientFluxes(
            "Cannot calculate quasi-bolometric flux with fewer " +
            "than two fluxes in the SED, {0} received".format(len(SED)))
    
    return QuasiBolometricFlux(
            value = integral_calculator.calculate(SED),
            uncertainty = uncertainty_calculator(SED),
            time = SED[0].time)

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

# TODO No test written
def calculate_qbol_flux(flux_group):
    """Turn a group of fluxes into a quasi-bolometric flux"""
    return get_quasi_bolometric_flux(TrapezoidalIntegralCalculator(),
                                     uncertainty_calculator_trapezoidal,
                                     flux_group)
