import unittest
import math
import numpy as np

from unittest.mock import Mock

from .context import superbol
from superbol import mag2flux
from superbol import fqbol

class TestGetQuasiBolometricFlux(unittest.TestCase):
    def setUp(self):
        self.integral_calculator = Mock()
        self.uncertainty_calculator = Mock()
        self.time = 1234.5

    def test_no_fluxes(self):
        with self.assertRaises(fqbol.InsufficientFluxes):
            fqbol.get_quasi_bolometric_flux(
                integral_calculator = self.integral_calculator,
                uncertainty_calculator = self.uncertainty_calculator,
                SED=[])

    def test_one_flux(self):
        flux = mag2flux.MonochromaticFlux(flux = 200,
                                          flux_uncertainty = 30,
                                          wavelength = 1,
                                          time = self.time)

        with self.assertRaises(fqbol.InsufficientFluxes):
            fqbol.get_quasi_bolometric_flux(
                integral_calculator = self.integral_calculator,
                uncertainty_calculator = self.uncertainty_calculator,
                SED=[flux])

    def test_two_fluxes(self):
        flux1 = mag2flux.MonochromaticFlux(flux = 100,
                                           flux_uncertainty = 0,
                                           wavelength = 0,
                                           time = self.time)

        flux2 = mag2flux.MonochromaticFlux(flux = 100,
                                           flux_uncertainty = 0,
                                           wavelength= 1,
                                           time = self.time)

        two_fluxes = [flux1, flux2]

        self.integral_calculator.calculate = Mock(return_value = 100)
        self.uncertainty_calculator = Mock(return_value = 10)
        
        expected_value = 100
        expected_uncertainty = 10

        result = fqbol.get_quasi_bolometric_flux(
            integral_calculator = self.integral_calculator,
            uncertainty_calculator = self.uncertainty_calculator,
            SED = two_fluxes)

        self.uncertainty_calculator.assert_called_once_with(two_fluxes)
        self.integral_calculator.calculate.assert_called_once_with(two_fluxes)
        self.assertEqual(expected_value, result.value)
        self.assertEqual(expected_uncertainty, result.uncertainty)
        self.assertEqual(self.time, result.time)

class TestUncertaintyCalculatorSpline(unittest.TestCase):

    def test_equal_flux_zero_uncertainty(self):
        flux1 = mag2flux.MonochromaticFlux(flux = 100,
                                           flux_uncertainty = 0,
                                           wavelength = 0,
                                           time = 0)

        flux2 = mag2flux.MonochromaticFlux(flux = 100,
                                           flux_uncertainty = 0,
                                           wavelength= 1,
                                           time = 0)

        expected_uncertainty = 0
        result = fqbol.uncertainty_calculator_spline(fluxes=[flux1, flux2])
        self.assertEqual(expected_uncertainty, result)

    def test_equal_flux_equal_uncertainty(self):
        flux1 = mag2flux.MonochromaticFlux(flux = 100,
                                           flux_uncertainty = 10,
                                           wavelength = 0,
                                           time = 0)

        flux2 = mag2flux.MonochromaticFlux(flux = 100,
                                           flux_uncertainty = 10,
                                           wavelength= 1,
                                           time = 0)

        expected_uncertainty = 10
        result = fqbol.uncertainty_calculator_spline(fluxes=[flux1,flux2])
        self.assertEqual(expected_uncertainty, result)

    def test_unequal_flux_unequal_uncertainty(self):
        flux1 = mag2flux.MonochromaticFlux(flux = 100,
                                           flux_uncertainty = 10,
                                           wavelength = 0,
                                           time = 0)

        flux2 = mag2flux.MonochromaticFlux(flux = 200,
                                           flux_uncertainty = 20,
                                           wavelength= 1,
                                           time = 0)

        expected_uncertainty = 15
        result = fqbol.uncertainty_calculator_spline(fluxes=[flux1,flux2])
        self.assertEqual(expected_uncertainty, result)

    def test_three_fluxes(self):
        flux1 = mag2flux.MonochromaticFlux(flux = 100,
                                           flux_uncertainty = 10,
                                           wavelength = 0,
                                           time = 0)

        flux2 = mag2flux.MonochromaticFlux(flux = 200,
                                           flux_uncertainty = 20,
                                           wavelength= 1,
                                           time = 0)

        flux3 = mag2flux.MonochromaticFlux(flux = 150,
                                           flux_uncertainty = 8,
                                           wavelength= 2,
                                           time = 0)

        expected_uncertainty = 31.75
        result = fqbol.uncertainty_calculator_spline(fluxes=[flux1,flux2,flux3])
        self.assertEqual(expected_uncertainty, result)

class TestSplineIntegralCalculator(unittest.TestCase):
    
    def setUp(self):
        self.integral_calculator = fqbol.SplineIntegralCalculator()
        self.flux1 = mag2flux.MonochromaticFlux(flux = 0,
                                                flux_uncertainty = 10,
                                                wavelength = 0,
                                                time = 0)
        self.flux2 = mag2flux.MonochromaticFlux(flux = 0.5,
                                                flux_uncertainty = 20,
                                                wavelength = 1,
                                                time = 0)
        self.flux3 = mag2flux.MonochromaticFlux(flux = 2,
                                                flux_uncertainty = 8,
                                                wavelength = 2,
                                                time = 0)
        self.flux4 = mag2flux.MonochromaticFlux(flux = 1.5,
                                                flux_uncertainty = 0.1,
                                                wavelength = 3,
                                                time = 0)
        self.fluxes = [self.flux1, self.flux2, self.flux3, self.flux4]

    def test_sort_fluxes(self):
        fluxes = [self.flux2, self.flux1, self.flux3, self.flux4]
        expected = [self.flux1, self.flux2, self.flux3, self.flux4]
        self.integral_calculator._sort_fluxes_by_wavelength(fluxes) 
        self.assertEqual(fluxes, expected)
    
    def test_flux_list(self):
        flux_list = self.integral_calculator._get_flux_list(self.fluxes)
        self.assertEqual([0, 0.5, 2, 1.5], flux_list)

    def test_wavelength_list(self):
        wl_list = self.integral_calculator._get_wavelength_list(self.fluxes)
        self.assertEqual([0, 1, 2, 3], wl_list)

    def test_spline_integral(self):
        integral = self.integral_calculator.calculate(self.fluxes)
        self.assertAlmostEqual(3.35, integral)
