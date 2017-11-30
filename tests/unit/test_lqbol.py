import unittest
import math
import numpy as np

from unittest.mock import Mock

from .context import superbol
from superbol import mag2flux
from superbol import lqbol

class TestGetQuasiBolometricFlux(unittest.TestCase):
    def setUp(self):
        self.integral_calculator = Mock()
        self.uncertainty_calculator = Mock()
        self.time = 1234.5

    def test_no_fluxes(self):
        with self.assertRaises(lqbol.InsufficientFluxes):
            lqbol.get_quasi_bolometric_flux(
                integral_calculator = self.integral_calculator,
                uncertainty_calculator = self.uncertainty_calculator,
                fluxes=[])

    def test_one_flux(self):
        flux = mag2flux.MonochromaticFlux(flux = 200,
                                          flux_uncertainty = 30,
                                          wavelength = 1,
                                          time = self.time)

        with self.assertRaises(lqbol.InsufficientFluxes):
            lqbol.get_quasi_bolometric_flux(
                integral_calculator = self.integral_calculator,
                uncertainty_calculator = self.uncertainty_calculator,
                fluxes=[flux])

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

        result = lqbol.get_quasi_bolometric_flux(
            integral_calculator = self.integral_calculator,
            uncertainty_calculator = self.uncertainty_calculator,
            fluxes = two_fluxes)

        self.uncertainty_calculator.assert_called_once_with(two_fluxes)
        self.integral_calculator.calculate.assert_called_once_with(two_fluxes)
        self.assertEqual(expected_value, result.value)
        self.assertEqual(expected_uncertainty, result.uncertainty)
        self.assertEqual(self.time, result.time)

class TestUncertaintyCalculatorTrapezoidal(unittest.TestCase):

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
        result = lqbol.uncertainty_calculator_trapezoidal(fluxes=[flux1, flux2])
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

        expected_uncertainty = math.sqrt(50)
        result = lqbol.uncertainty_calculator_trapezoidal(fluxes=[flux1,flux2])
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

        expected_uncertainty = math.sqrt(125)
        result = lqbol.uncertainty_calculator_trapezoidal(fluxes=[flux1,flux2])
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

        expected_uncertainty = math.sqrt(441)
        result = lqbol.uncertainty_calculator_trapezoidal(fluxes=[flux1,flux2,flux3])
        self.assertEqual(expected_uncertainty, result)

class TestTrapezoidalIntegralCalculator(unittest.TestCase):
    
    def setUp(self):
        self.integral_calculator = lqbol.TrapezoidalIntegralCalculator()
        self.flux1 = mag2flux.MonochromaticFlux(flux = 100,
                                                flux_uncertainty = 10,
                                                wavelength = 0,
                                                time = 0)
        self.flux2 = mag2flux.MonochromaticFlux(flux = 200,
                                                flux_uncertainty = 20,
                                                wavelength= 1,
                                                time = 0)
        self.flux3 = mag2flux.MonochromaticFlux(flux = 150,
                                                flux_uncertainty = 8,
                                                wavelength= 2,
                                                time = 0)
        self.fluxes = [self.flux1, self.flux2, self.flux3]

    def test_sort_fluxes(self):
        fluxes = [self.flux2, self.flux1, self.flux3]
        expected = [self.flux1, self.flux2, self.flux3]
        self.integral_calculator._sort_fluxes_by_wavelength(fluxes) 
        self.assertEqual(fluxes, expected)
    
    def test_flux_list(self):
        flux_list = self.integral_calculator._get_flux_list(self.fluxes)
        self.assertEqual([100, 200, 150], flux_list)

    def test_wavelength_list(self):
        wl_list = self.integral_calculator._get_wavelength_list(self.fluxes)
        self.assertEqual([0, 1, 2], wl_list)

    def test_trapezoidal_integral(self):
        integral = self.integral_calculator.calculate(self.fluxes)
        self.assertEqual(325, integral)

class TestFluxLuminosityConverter(unittest.TestCase):

    def setUp(self):
        self.fqbol = lqbol.QuasiBolometricFlux(value = 10, uncertainty = 1, time = 0)
        self.distance = lqbol.Distance(value = 100, uncertainty = 10)

    def test_convert_flux_to_luminosity(self):
        expected = self.fqbol.value * 4.0 * math.pi * self.distance.value**2
        result = lqbol.convert_flux_to_luminosity(self.fqbol, self.distance)
        self.assertEqual(expected, result.value)

    def test_convert_flux_to_luminosity_uncertainty(self):
        expected = math.sqrt((4.0 * math.pi * self.distance.value**2 * self.fqbol.uncertainty)**2 + (self.fqbol.value * 8.0 * math.pi * self.distance.value * self.distance.uncertainty)**2)
        result = lqbol.convert_flux_to_luminosity(self.fqbol, self.distance)
        self.assertEqual(expected, result.uncertainty)

    def test_convert_flux_to_luminosity_time(self):
        expected = 0
        result = lqbol.convert_flux_to_luminosity(self.fqbol, self.distance)
        self.assertEqual(expected, result.time)
