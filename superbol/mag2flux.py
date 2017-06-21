import math

class ObservedMagnitude(object):

    def __init__(self, magnitude, uncertainty, band):
        self.magnitude = magnitude
        self.uncertainty = uncertainty
        self.band = band

    def convert_to_flux(self):
        monochromatic_flux = MagnitudeToFluxConverter(self).convert_to_flux()
        return monochromatic_flux

class Band(object):

    def __init__(self, name, effective_wavelength, flux_conversion_factor):
        self.name = name
        self.effective_wavelength = effective_wavelength
        self.flux_conversion_factor = flux_conversion_factor

class MonochromaticFlux(object):

    def __init__(self, flux, flux_uncertainty, wavelength):
        self.flux = flux
        self.flux_uncertainty = flux_uncertainty
        self.wavelength = wavelength

class MagnitudeToFluxConverter(object):

    def __init__(self, observation):
        self.observation = observation

    def _calculate_flux(self):
        flux_conversion_factor = self.observation.band.flux_conversion_factor
        flux = flux_conversion_factor * 10**(-0.4 * self.observation.magnitude)
        return flux

    def _calculate_flux_uncertainty(self, flux):
        flux_uncertainty = flux * 0.4 * math.log(10) * self.observation.uncertainty
        return flux_uncertainty

    def convert_to_flux(self):
        flux = self._calculate_flux()
        flux_uncertainty = self._calculate_flux_uncertainty(flux)
        effective_wavelength = self.observation.band.effective_wavelength
        monochromatic_flux = MonochromaticFlux(flux, flux_uncertainty, effective_wavelength)
        return monochromatic_flux
