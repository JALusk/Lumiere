import unittest

from .context import superbol
from superbol import faug
from superbol import mag2flux
from superbol import blackbody

class TestTrimSED(unittest.TestCase):

    def setUp(self):
        self.time = 1234.5
        self.flux1 = mag2flux.MonochromaticFlux(flux = 1,
                                                flux_uncertainty = 1,
                                                wavelength = 3660.0,
                                                time = self.time)
        self.flux2 = mag2flux.MonochromaticFlux(flux = 2,
                                                flux_uncertainty = 1,
                                                wavelength = 4380.0,
                                                time = self.time)
        self.flux3 = mag2flux.MonochromaticFlux(flux = 8,
                                                flux_uncertainty = 1,
                                                wavelength = 5450.0,
                                                time = self.time)
        self.flux4 = mag2flux.MonochromaticFlux(flux = 4,
                                                flux_uncertainty = 1,
                                                wavelength = 6410.0,
                                                time = self.time)
        self.flux5 = mag2flux.MonochromaticFlux(flux = 5,
                                                flux_uncertainty = 1,
                                                wavelength = 7980.0,
                                                time = self.time)
        self.SED = [self.flux1, self.flux2, self.flux3, self.flux4, self.flux5]

    def test_trim_SED_5000_angstroms(self):
        expected = [self.flux3, self.flux4, self.flux5]
        result = faug.trim_SED(self.SED, 5000.0)
        
        self.assertEqual(expected, result)

    def test_trim_SED_7000_angstroms(self):
        expected = [self.flux5]
        result = faug.trim_SED(self.SED, 7000.0)

        self.assertEqual(expected, result)

    # Finds the max flux from SED
    def test_find_max_flux(self):
        self.flux4.flux = 10    # temporary spike in flux to test the max function
        expected = self.flux4 # change to max(self.SED)
        result = faug.find_max_flux(self.SED)
        self.assertEqual(expected, result)

    # Gets the trimmed SED; all fluxes with wavelength >= the max flux- (3, 4, 5)
    def test_get_max_flux_wavelength(self):
        ''' Test if trimming left of max flux works'''
        max_flux = faug.find_max_flux(self.SED) # gets max flux
        max_flux_wavelength = max_flux.wavelength   # max flux's wavelength

        expected = [self.flux3, self.flux4, self.flux5] # This is what should remain after trim
        result = faug.trim_SED(self.SED, max_flux_wavelength)   
        self.assertEqual(expected, result)

    def test_equal_max_fluxes(self):
        ''' If there are 2 'max' fluxes, return the first one'''
        self.flux3.flux = 10    # Set flux 3 and 4 to be max flux
        self.flux4.flux = 10
        expected = self.flux3   # Should return the first in the list (sorted by wavelength)
        result = faug.find_max_flux(self.SED)
        self.assertEqual(expected, result)

    def test_trim_SED_to_peak(self):
        ''' Keeps only the fluxes in SED that have a wavelength >= max_flux.wavelength'''
        expected = [self.flux3, self.flux4, self.flux5]
        result = faug.trim_SED_to_peak(self.SED)
        self.assertEqual(expected, result) 
    
    def test_find_min_flux(self):
        self.flux1.flux = 1    # temporary spike in flux to test the min function
        expected = self.flux1
        result = faug.find_min_flux(self.SED)
        self.assertEqual(expected, result)


# TODO Add tests

#class TestIRCorrection(unittest.TestCase):
#    pass


