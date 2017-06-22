import math

class Band(object):
    
    def __init__(self, name, effective_wavelength, flux_conversion_factor):
        self.name = name
        self.effective_wavelength = effective_wavelength
        self.flux_conversion_factor = flux_conversion_factor

class ObservedMagnitude(object):

    def __init__(self, magnitude, uncertainty, band):
        self.magnitude = magnitude
        self.uncertainty = uncertainty
        self.band = band

    def convert_to_flux(self):
        monochromatic_flux = MagnitudeToFluxConverter(self).convert()
        return monochromatic_flux

class MonochromaticFlux(object):

    def __init__(self, flux, flux_uncertainty, wavelength):
        self.flux = flux
        self.flux_uncertainty = flux_uncertainty
        self.wavelength = wavelength

class MagnitudeToFluxConverter(object):

    def __init__(self, observed_magnitude):
        self.observed_magnitude = observed_magnitude

    def _calculate_flux(self):
        flux_conversion_factor = self.observed_magnitude.band.flux_conversion_factor
        flux = flux_conversion_factor * 10**(-0.4 * self.observed_magnitude.magnitude)
        return flux

    def _calculate_flux_uncertainty(self, flux):
        flux_uncertainty = flux * 0.4 * math.log(10) * self.observed_magnitude.uncertainty
        return flux_uncertainty

    def convert(self):
        flux = self._calculate_flux()
        flux_uncertainty = self._calculate_flux_uncertainty(flux)
        effective_wavelength = self.observed_magnitude.band.effective_wavelength
        monochromatic_flux = MonochromaticFlux(flux, flux_uncertainty, effective_wavelength)
        return monochromatic_flux
