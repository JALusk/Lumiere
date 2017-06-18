import math

class Observation(object):

    def __init__(self, magnitude, uncertainty, flux_conversion_factor):
        self.magnitude = magnitude
        self.uncertainty = uncertainty
        self.flux_conversion_factor = flux_conversion_factor

    def convert_to_flux(self):
        flux = self.flux_conversion_factor * 10**(-0.4 * self.magnitude)
        return flux

    def calculate_flux_uncertainty(self):
        flux_uncertainty = self.convert_to_flux() * 0.4 * math.log(10) * self.uncertainty
        return flux_uncertainty
