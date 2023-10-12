import unittest

from .context import superbol
from superbol import faug
from superbol import mag2flux
from superbol.blackbody import *
from superbol import fqbol

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

    def test_trim_SED_no_wavelength(self):
        expected = [self.flux3, self.flux4, self.flux5]
        result = faug.trim_SED(self.SED)
        
        self.assertEqual(expected, result)

    def test_find_max_flux(self):
        self.flux4.flux = 10  
        expected = self.flux4
        result = faug.find_max_flux(self.SED)
        self.assertEqual(expected, result)

    def test_find_longest_wavelength_is_flux5(self):
        self.flux5.wavelength = 10000
        expected = self.flux5.wavelength
        result = faug.find_longest_wavelength(self.SED)
        self.assertEqual(expected, result)

    def test_find_longest_wavelength_is_flux3(self):
        self.flux3.wavelength = 10000
        expected = self.flux3.wavelength
        result = faug.find_longest_wavelength(self.SED)
        self.assertEqual(expected, result)

    def test_find_longest_wavelength_is_flux1(self):
        self.flux1.wavelength = 10000
        expected = self.flux1.wavelength
        result = faug.find_longest_wavelength(self.SED)
        self.assertEqual(expected, result)

    def test_find_shortest_wavelength_is_flux5(self):
        self.flux5.wavelength = 1
        expected = self.flux5.wavelength
        result = faug.find_shortest_wavelength(self.SED)
        self.assertEqual(expected, result)

    def test_find_shortest_wavelength_is_flux3(self):
        self.flux3.wavelength = 1
        expected = self.flux3.wavelength
        result = faug.find_shortest_wavelength(self.SED)
        self.assertEqual(expected, result)

    def test_find_shortest_wavelength_is_flux1(self):
        self.flux1.wavelength = 1
        expected = self.flux1.wavelength
        result = faug.find_shortest_wavelength(self.SED)
        self.assertEqual(expected, result)

    def test_get_max_flux_wavelength(self):
        ''' Test if trimming left of max flux works'''
        max_flux = faug.find_max_flux(self.SED)
        max_flux_wavelength = max_flux.wavelength   

        expected = [self.flux3, self.flux4, self.flux5]
        result = faug.trim_SED(self.SED, max_flux_wavelength)   
        self.assertEqual(expected, result)

    def test_equal_max_fluxes(self):
        ''' If there are 2 'max' fluxes, return the first one'''
        self.flux3.flux = 10    
        self.flux4.flux = 10
        expected = self.flux3
        result = faug.find_max_flux(self.SED)
        self.assertEqual(expected, result)

    def test_trim_SED_to_peak(self):
        ''' Keeps only the fluxes in SED that have a wavelength >= max_flux.wavelength'''
        expected = [self.flux3, self.flux4, self.flux5]
        result = faug.trim_SED_to_peak(self.SED)
        self.assertEqual(expected, result) 
    
    def test_find_bluest_flux(self):
        self.flux1.flux = 1
        expected = self.flux1
        result = faug.find_bluest_flux(self.SED)
        self.assertEqual(expected, result)

    def test_get_IR(self):
        self.flux5.wavelength = 10000
        bbfit = BlackbodyFit()
        bbfit.fit_to_SED(self.SED)
        expected = (bb_total_flux(bbfit.temperature, bbfit.angular_radius) 
                    - bb_flux_integrated(self.flux5.wavelength, bbfit.temperature, bbfit.angular_radius))
        result = faug.get_IR(self.SED)

        self.assertEqual(expected, result)

class TestTrimRealSED(unittest.TestCase):
    
    def setUp(self):
        # Read data from SN 2018hna
        self.time = 2458493.0
        self.flux1 = mag2flux.MonochromaticFlux(flux = 5.354E-15,
                                                flux_uncertainty = 1.479E-16,
                                                wavelength = 4380.0,
                                                time = self.time)
        self.flux2 = mag2flux.MonochromaticFlux(flux = 7.586E-15,
                                                flux_uncertainty = 3.494E-16,
                                                wavelength = 5450.0,
                                                time = self.time)
        self.flux3 = mag2flux.MonochromaticFlux(flux = 6.593E-15,
                                                flux_uncertainty = 1.214E-16,
                                                wavelength = 6410.0,
                                                time = self.time)
        self.flux4 = mag2flux.MonochromaticFlux(flux = 4.040E-15,
                                                flux_uncertainty = 1.488E-16,
                                                wavelength = 7980.0,
                                                time = self.time)
        self.flux5 = mag2flux.MonochromaticFlux(flux = 1.503E-15,
                                                flux_uncertainty = 4.153E-17,
                                                wavelength = 12200.0,
                                                time = self.time)
        self.flux6 = mag2flux.MonochromaticFlux(flux = 6.701E-16,
                                                flux_uncertainty = 1.234E-17,
                                                wavelength = 16300.0,
                                                time = self.time)
        self.SED = [self.flux1, self.flux2, self.flux3, self.flux4, self.flux5, self.flux6]

    def test_trim_real_SED_5000_angstroms(self):
        expected = [self.flux2, self.flux3, self.flux4, self.flux5, self.flux6]
        result = faug.trim_SED(self.SED, 5000.0)
        
        self.assertEqual(expected, result)

    def test_trim_real_SED_7000_angstroms(self):
        expected = [self.flux4, self.flux5, self.flux6]
        result = faug.trim_SED(self.SED, 7000.0)
        
        self.assertEqual(expected, result)

    def test_find_real_SED_max_flux(self):
        expected = self.flux2
        result = faug.find_max_flux(self.SED)

        self.assertEqual(expected, result)

    def test_trim_real_SED_to_peak(self):
        expected = [self.flux2, self.flux3, self.flux4, self.flux5, self.flux6]
        result = faug.trim_SED_to_peak(self.SED)
        self.assertEqual(expected, result)

    def test_find_real_SED_min_flux(self):
        expected = self.flux1
        result = faug.find_bluest_flux(self.SED)
        self.assertEqual(expected, result)


class Test_UV_IR_Corrections(unittest.TestCase):
    def setUp(self):
        # Read data from SN 2018hna
        self.time = 2458493.0
        self.flux1 = mag2flux.MonochromaticFlux(flux = 5.354E-15,
                                                flux_uncertainty = 1.479E-16,
                                                wavelength = 4380.0,
                                                time = self.time)
        self.flux2 = mag2flux.MonochromaticFlux(flux = 7.586E-15,
                                                flux_uncertainty = 3.494E-16,
                                                wavelength = 5450.0,
                                                time = self.time)
        self.flux3 = mag2flux.MonochromaticFlux(flux = 6.593E-15,
                                                flux_uncertainty = 1.214E-16,
                                                wavelength = 6410.0,
                                                time = self.time)
        self.flux4 = mag2flux.MonochromaticFlux(flux = 4.040E-15,
                                                flux_uncertainty = 1.488E-16,
                                                wavelength = 7980.0,
                                                time = self.time)
        self.flux5 = mag2flux.MonochromaticFlux(flux = 1.503E-15,
                                                flux_uncertainty = 4.153E-17,
                                                wavelength = 12200.0,
                                                time = self.time)
        self.flux6 = mag2flux.MonochromaticFlux(flux = 6.701E-16,
                                                flux_uncertainty = 1.234E-17,
                                                wavelength = 16300.0,
                                                time = self.time)
        self.SED = [self.flux1, self.flux2, self.flux3, self.flux4, self.flux5, self.flux6]

    def test_fit_blackbody(self):
        test_bbfit = BlackbodyFit()
        test_trimSED = faug.trim_SED(self.SED) 
        test_bbfit.fit_to_SED(test_trimSED)

        self.assertAlmostEqual(test_bbfit.temperature, 6421.8, 1)

    def test_get_UV(self):
        base = (self.flux1.wavelength - 2000)
        height = self.flux1.flux
        expected = (0.5 * (base * height))
        result = faug.get_UV(self.SED)

        self.assertAlmostEqual(expected, result)
        self.assertEqual(expected, result)
