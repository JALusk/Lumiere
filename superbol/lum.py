import math

def convert_flux_to_luminosity(flux, distance):
    """Convert bolometric flux to bolometric luminosity"""
    lbol = 4.0 * math.pi * distance.value**2 * flux.value
    uncertainty = math.sqrt((2.0 * lbol / distance.value * distance.uncertainty)**2
                            + (4.0 * math.pi * distance.value**2 * flux.uncertainty)**2)
    return BolometricLuminosity(lbol, uncertainty)

class BolometricFlux(object):

    def __init__(self, value, uncertainty):
        if value < 0:
            raise ValueError("Negative luminosity")
        else:
            self.value = value
            self.uncertainty = abs(uncertainty)

    def to_lbol(self, distance):
        """Convert bolometric flux to bolometric luminosity"""
        return convert_flux_to_luminosity(self, distance)

class Distance(object):

    def __init__(self, value, uncertainty):
        if value < 0:
            raise ValueError("Negative distance")
        else:
            self.value = value
            self.uncertainty = abs(uncertainty)

class BolometricLuminosity(object):

    def __init__(self, value, uncertainty):
        if value < 0:
            raise ValueError("Negative luminosity")
        else:
            self.value = value
            self.uncertainty = abs(uncertainty)
