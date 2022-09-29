import math
import numpy as np
from matplotlib import pyplot as plt
import scipy

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

class SplineIntegralCalculator(object):
    """Integrate between fluxes using a cubic spline"""
    def calculate(self, fluxes):
        """Take a cubic spline of the fluxes and wavelengths and integrate"""
        self._sort_fluxes_by_wavelength(fluxes)
        flux_list = self._get_flux_list(fluxes)
        wavelength_list = self._get_wavelength_list(fluxes)
        print("Wavelength list: ", wavelength_list)
        print("Flux list: ", flux_list)
        spline = scipy.interpolate.CubicSpline(wavelength_list, flux_list, bc_type='natural') 
        #print([x for x in spline])
        integrated_spline = float(spline.integrate(wavelength_list[0], wavelength_list[-1]))
        
        #Plots spline for reference or poster
        new_wavel = np.linspace(wavelength_list[0], wavelength_list[-1], 200)
        spline_flux = spline(new_wavel)
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("u(wavelength) (kJ/nm)")
        plt.scatter(wavelength_list, flux_list)
        plt.plot(new_wavel, spline_flux)
        plt.savefig('../Desktop/spline_qbol.png')
        #nope gotta plot magnitude not flux

        return integrated_spline

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


def uncertainty_calculator_spline(fluxes):
    """Take cubic spline of the flux uncertainties"""
    
    def _get_flux_uncertainty_list(fluxes):
        """Return a list of the flux uncertainties"""
        return [f.flux_uncertainty for f in fluxes]
    
    def _get_wavelength_list(fluxes):
        """Return a list of flux wavelengths"""
        return [f.wavelength for f in fluxes]

    def _get_flux_list(fluxes):
        """Return a list of flux wavelengths"""
        return [f.flux for f in fluxes]
    
    wavelength_list = _get_wavelength_list(fluxes)
    flux_uncertainty_list = _get_flux_uncertainty_list(fluxes)
    flux_list = _get_flux_list(fluxes)
    flux_plus_uncertainty = np.add(flux_list, flux_uncertainty_list)

    uncertainty_spline = scipy.interpolate.CubicSpline(wavelength_list, flux_plus_uncertainty, bc_type='natural')
    uncertainty_integrated = float(uncertainty_spline.integrate(wavelength_list[0], wavelength_list[-1]))

    flux_spline = scipy.interpolate.CubicSpline(wavelength_list, flux_list, bc_type='natural') 
    flux_integrated = float(flux_spline.integrate(wavelength_list[0], wavelength_list[-1]))
    
    qbolflux_uncertainty = uncertainty_integrated - flux_integrated
    #ratio_uncertainty = qbolflux_uncertainty / flux_integrated

    return qbolflux_uncertainty

def calculate_qbol_flux(flux_group):
    """Turn a group of fluxes into a quasi-bolometric flux"""
    print(flux_group[0].time)
    return get_quasi_bolometric_flux(SplineIntegralCalculator(),
                                     uncertainty_calculator_spline,
                                     flux_group)