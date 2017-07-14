import math

class Band(object):
    
    def __init__(self, name, effective_wavelength, flux_conversion_factor):
        self.name = name
        self.effective_wavelength = effective_wavelength
        self.flux_conversion_factor = flux_conversion_factor

class ObservedMagnitude(object):

    def __init__(self, magnitude, uncertainty, band, time):
        self.magnitude = magnitude
        self.uncertainty = uncertainty
        self.band = band
        self.time = time

    def convert_to_flux(self):
        monochromatic_flux = MagnitudeToFluxConverter().convert(self)
        return monochromatic_flux

class MonochromaticFlux(object):

    def __init__(self, flux, flux_uncertainty, wavelength, time):
        self.flux = flux
        self.flux_uncertainty = flux_uncertainty
        self.wavelength = wavelength
        self.time = time

class MagnitudeToFluxConverter(object):

    def _calculate_flux(self, magnitude, flux_conversion_factor):
        flux = flux_conversion_factor * 10**(-0.4 * magnitude)
        return flux

    def _calculate_flux_uncertainty(self, flux, magnitude_uncertainty):
        flux_uncertainty = flux * 0.4 * math.log(10) * magnitude_uncertainty
        return flux_uncertainty

    def convert(self, observed_magnitude):
        """Convert an observed magnitude to a monochromatic flux"""
        magnitude = observed_magnitude.magnitude
        uncertainty = observed_magnitude.uncertainty
        time = observed_magnitude.time
        flux_conversion_factor = observed_magnitude.band.flux_conversion_factor

        flux = self._calculate_flux(magnitude, flux_conversion_factor)
        flux_uncertainty = self._calculate_flux_uncertainty(flux, uncertainty)
        wavelength = observed_magnitude.band.effective_wavelength
        monochromatic_flux = MonochromaticFlux(flux, flux_uncertainty, wavelength, time)
        return monochromatic_flux
