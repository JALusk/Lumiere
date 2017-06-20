import math

class Observation(object):

    def __init__(self, magnitude, uncertainty, band):
        self.magnitude = magnitude
        self.uncertainty = uncertainty
        self.band = band

    def convert_to_flux(self):
        flux = self.band.flux_conversion_factor * 10**(-0.4 * self.magnitude)
        flux_uncertainty = flux * 0.4 * math.log(10) * self.uncertainty
        wavelength = self.band.effective_wavelength
        monochromatic_flux = MonochromaticFlux(flux, flux_uncertainty, wavelength)
        return monochromatic_flux

class Band(object):

    def __init__(self, name, effective_wavelength, flux_conversion_factor):
        self.name = name
        self.effective_wavelength = effective_wavelength
        self.flux_conversion_factor = flux_conversion_factor

class MonochromaticFlux(object):

    def __init__(self, flux, uncertainty, wavelength):
        self.flux = flux
        self.uncertainty = uncertainty
        self.wavelength = wavelength
