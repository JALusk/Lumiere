class Observation(object):

    def __init__(self, magnitude):
        self.magnitude = magnitude

    def convert_to_flux(self):
        return 417.5E-11
