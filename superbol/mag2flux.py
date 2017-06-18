class Observation(object):

    def __init__(self, magnitude, flux_conversion_factor):
        self.magnitude = magnitude
        self.flux_conversion_factor = flux_conversion_factor

    def convert_to_flux(self):
        flux = self.flux_conversion_factor * 10**(-0.4 * self.magnitude)
        return flux
