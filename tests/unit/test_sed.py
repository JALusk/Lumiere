import unittest
import numpy as np

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

    def test_get_integrable_fluxes(self):
        result = sed.get_integrable_fluxes(self.fluxes)
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
