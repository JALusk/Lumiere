import unittest
import numpy as np
import math

from .context import superbol
from superbol import sed
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
        result = sed.group_fluxes(fluxes)
        self.assertEqual(expected, result)

    def test_group_fluxes_floor_different_days(self):
        fluxes = [self.flux10, self.flux11, self.flux12, self.flux21, self.flux27]
        expected = [[self.flux10, self.flux11, self.flux12], 
                    [self.flux21, self.flux27]]
        result = sed.group_fluxes(fluxes)
        self.assertEqual(expected, result)

class TestCombineFluxes(unittest.TestCase):

    def setUp(self):
        self.flux1 = mag2flux.MonochromaticFlux(flux = 100,
                                                flux_uncertainty = 10,
                                                wavelength = 1,
                                                time = 0)
        self.flux2 = mag2flux.MonochromaticFlux(flux = 200,
                                                flux_uncertainty = 10,
                                                wavelength= 1,
                                                time = 0)
        self.flux3 = mag2flux.MonochromaticFlux(flux = 150,
                                                flux_uncertainty = 8,
                                                wavelength= 2,
                                                time = 0)
        self.flux4 = mag2flux.MonochromaticFlux(flux = 50,
                                                flux_uncertainty = 8,
                                                wavelength= 3,
                                                time = 0)
        self.flux5 = mag2flux.MonochromaticFlux(flux = 60,
                                                flux_uncertainty = 8,
                                                wavelength= 3,
                                                time = 0)
        self.fluxes = [self.flux1, self.flux2, self.flux3, self.flux4, self.flux5]
        self.repeated_fluxes1 = [self.flux1, self.flux2]
        self.repeated_fluxes3 = [self.flux4, self.flux5]

    def test_combine_fluxes_equal_uncertainties(self):
        result = sed.combine_fluxes(self.repeated_fluxes1)
        expected = mag2flux.MonochromaticFlux(flux = 150,
                                              flux_uncertainty = np.sqrt(200)/2.,
                                              wavelength = 1,
                                              time = 0)
        self.assertEqual(expected, result)

    def test_combine_fluxes_unequal_uncertainties(self):
        fluxes = [self.flux1, self.flux3]
        result = sed.combine_fluxes(fluxes)
        expected = mag2flux.MonochromaticFlux(flux = 130.488,
                                              flux_uncertainty = 6.247,
                                              wavelength = 1,
                                              time = 0)
        self.assertAlmostEqual(expected.flux, result.flux, 3)
        self.assertAlmostEqual(expected.flux_uncertainty, result.flux_uncertainty, 3)


    def test_yield_fluxes_at_each_observed_wavelength(self):
        result_generator = sed.yield_fluxes_at_each_observed_wavelength(self.fluxes)
        expected1 = self.repeated_fluxes1
        expected3 = self.repeated_fluxes3
        self.assertEqual(expected1, next(result_generator))
        self.assertEqual([self.flux3], next(result_generator))
        self.assertEqual(expected3, next(result_generator))
        with self.assertRaises(StopIteration):
            next(result_generator)

    def test_get_SED(self):
        result = sed.get_SED(self.fluxes)
        expected = [mag2flux.MonochromaticFlux(150, np.sqrt(200)/2., 1, 0), self.flux3,
                    mag2flux.MonochromaticFlux(55.0, np.sqrt(128)/2., 3, 0)]

        # Ugly
        for i, flux in enumerate(result):
            self.assertAlmostEqual(result[i].flux, expected[i].flux)
            self.assertAlmostEqual(result[i].flux_uncertainty, expected[i].flux_uncertainty)

    def test_weighted_average(self):
        expected = 10.4
        result = sed.weighted_average([10.0, 12.0], [0.5, 1.0])
        self.assertEqual(expected, result)

    def test_weighted_average_uncertainty(self):
        expected = 0.69
        uncertainties = [1, 1, 3]
        result = sed.weighted_average_uncertainty(uncertainties)
        self.assertAlmostEqual(expected, result, 2)

    def test_get_weights(self):
        expected = [4, 1]
        uncertainties = [0.5, 1.0]
        result = sed.get_weights(uncertainties)
        self.assertEqual(expected, result)

class TestInterpolateSED(unittest.TestCase):

    def setUp(self):
        self.flux01 = mag2flux.MonochromaticFlux(100, 2, 1, 0)
        self.flux02 = mag2flux.MonochromaticFlux(200, 2, 2, 0)
        self.flux03 = mag2flux.MonochromaticFlux(150, 2, 3, 0)
        self.flux11 = mag2flux.MonochromaticFlux(100, 2, 1, 1)
        self.flux13 = mag2flux.MonochromaticFlux(150, 2, 3, 1)
        self.flux21 = mag2flux.MonochromaticFlux(100, 2, 1, 2)
        self.flux22 = mag2flux.MonochromaticFlux(200, 2, 2, 2)
        self.flux23 = mag2flux.MonochromaticFlux(150, 2, 3, 2)
        self.flux31 = mag2flux.MonochromaticFlux(150, 2, 1, 3)
        self.flux32 = mag2flux.MonochromaticFlux(150, 2, 2, 3)
        self.SED0 = sed.get_SED([self.flux01, self.flux02, self.flux03])
        self.SED1 = sed.get_SED([self.flux11, self.flux13])
        self.SED2 = sed.get_SED([self.flux21, self.flux22, self.flux23])
        self.SED3 = sed.get_SED([self.flux31, self.flux32])

    def test_simple_interpolation(self):
        sed.interpolate_missing_fluxes([self.SED0, self.SED1, self.SED2])
        previous_flux = self.flux02
        next_flux = self.flux22
        unobserved_time = 1
        weight1 = (2-1)/(2-0)
        weight2 = (1-0)/(2-0)
        uncertainty = math.sqrt(weight1**2 * previous_flux.flux_uncertainty**2 + weight2**2 + next_flux.flux_uncertainty**2)
        interpolated_flux = mag2flux.MonochromaticFlux(200, uncertainty, 2, 1)
        self.assertTrue(interpolated_flux in self.SED1)

    def test_get_interpolated_fluxes(self):
        lightcurve = [self.flux02, self.flux22, self.flux32]
        previous_flux = self.flux02
        next_flux = self.flux22
        unobserved_time = 1
        weight1 = (2-1)/(2-0)
        weight2 = (1-0)/(2-0)
        uncertainty = math.sqrt(weight1**2 * previous_flux.flux_uncertainty**2 + weight2**2 + next_flux.flux_uncertainty**2)
        expected = [mag2flux.MonochromaticFlux(200.0, uncertainty , 2, 1)]
        result = sed.get_interpolated_fluxes(lightcurve, [1])
        self.assertEqual(expected, result)

    def test_append_interpolated_fluxes_to_SEDs(self):
        flux = mag2flux.MonochromaticFlux(200.0, 2, 2, 1)
        SED = self.SED1
        sed.append_interpolated_fluxes_to_SEDs([flux], [SED])
        self.assertTrue(flux in SED)

    def test_get_previous_flux(self):
        monochromatic_lc = [self.flux02, self.flux22, self.flux32]
        unobserved_time = 1
        previous_flux = sed.get_previous_flux(monochromatic_lc, unobserved_time)
        self.assertEqual(self.flux02, previous_flux)

    def test_get_next_flux(self):
        monochromatic_lc = [self.flux02, self.flux22, self.flux32]
        unobserved_time = 1
        next_flux = sed.get_next_flux(monochromatic_lc, unobserved_time)
        self.assertEqual(self.flux22, next_flux)

    def test_get_gap_size(self):
        times = [0, 1, 3, 6]
        unobserved_time = 2
        max_delta = sed.get_gap_size(times, unobserved_time)
        self.assertEqual(2, max_delta)

    def test_get_gap_size(self):
        times = [0, 1, 3, 6]
        unobserved_time = 4
        max_delta = sed.get_gap_size(times, unobserved_time)
        self.assertEqual(3, max_delta)

    def test_get_gap_size_unsorted(self):
        times = [0, 6, 3, 1]
        unobserved_time = 4
        max_delta = sed.get_gap_size(times, unobserved_time)
        self.assertEqual(3, max_delta)

    def test_big_gap(self):
        monochromatic_lc = [self.flux01, self.flux31]
        interpolated_lc = sed.get_interpolated_fluxes(monochromatic_lc, [0, 1, 2, 3])
        self.assertEqual([], interpolated_lc)

    def test_get_interpolated_flux_uncertainty(self):
        previous_flux = self.flux02
        next_flux = self.flux22
        unobserved_time = 1
        weight1 = (2-1)/(2-0)
        weight2 = (1-0)/(2-0)
        expected = math.sqrt(weight1**2 * previous_flux.flux_uncertainty**2 + weight2**2 + next_flux.flux_uncertainty**2)
        result = sed.get_interpolated_flux_uncertainty(previous_flux, next_flux, unobserved_time)
        self.assertEqual(expected, result)

    def test_get_unobserved_times(self):
        lightcurve = [self.flux02, self.flux22, self.flux32]
        observed_times = sed.get_observed_times([self.SED0, self.SED1, self.SED2, self.SED3])
        unobserved_times = sed.get_unobserved_times(lightcurve, observed_times)
        self.assertEqual([1], unobserved_times)
        
    def test_get_observed_wavelengths(self):
        expected = [1, 2, 3]
        result = sed.get_observed_wavelengths([self.SED0, self.SED1, self.SED2])
        self.assertEqual(expected, result)

    def test_get_observed_times(self):
        expected = [0, 1, 2]
        result = sed.get_observed_times([self.SED0, self.SED1, self.SED2])
        self.assertEqual(expected, result)

    def test_get_monochromatic_lightcurve(self):
        expected = [self.flux01, self.flux11, self.flux21]
        result = sed.get_monochromatic_lightcurve([self.SED0, self.SED1, self.SED2], 1)
        self.assertEqual(expected, result)

