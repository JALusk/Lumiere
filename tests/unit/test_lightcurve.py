import unittest
import math

from .context import superbol
from superbol import lightcurve
from superbol import mag2flux

class TestGroupFluxes(unittest.TestCase):

    def setUp(self):
        self.flux10 = mag2flux.MonochromaticFlux(flux = 0,
                                                 flux_uncertainty = 0,
                                                 wavelength = 0,
                                                 time = 1)

        self.flux11 = mag2flux.MonochromaticFlux(flux = 0,
                                                 flux_uncertainty = 0,
                                                 wavelength= 0,
                                                 time = 1.1)

        self.flux12 = mag2flux.MonochromaticFlux(flux = 0,
                                                 flux_uncertainty = 0,
                                                 wavelength= 0,
                                                 time = 1.3)

        self.flux21 = mag2flux.MonochromaticFlux(flux = 0,
                                                 flux_uncertainty = 0,
                                                 wavelength= 0,
                                                 time = 2.1)
        
        self.flux27 = mag2flux.MonochromaticFlux(flux = 0,
                                                 flux_uncertainty = 0,
                                                 wavelength= 0,
                                                 time = 2.7)

    def test_group_fluxes_floor_same_day(self):
        fluxes = [self.flux10, self.flux11, self.flux12]
        expected = [[self.flux10, self.flux11, self.flux12]]
        result = lightcurve.group_fluxes(fluxes, math.floor)
        self.assertEqual(expected, result)

    def test_group_fluxes_floor_different_days(self):
        fluxes = [self.flux10, self.flux11, self.flux12, self.flux21, self.flux27]
        expected = [[self.flux10, self.flux11, self.flux12], 
                    [self.flux21, self.flux27]]
        result = lightcurve.group_fluxes(fluxes, math.floor)
        self.assertEqual(expected, result)
