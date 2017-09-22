import unittest
import math
from .context import superbol
from superbol import mag2flux

class TestObservedMagnitude(unittest.TestCase):

    def setUp(self):
        self.B_band = mag2flux.Band('B', 4380.0, 632.0E-11)
        self.B_obs = mag2flux.ObservedMagnitude(18.793, 0.02, self.B_band,
                                                2451663.30)
        
    def test_convert_to_flux_converts_Johnson_B_magnitude_to_flux(self):
        expected_flux = self.B_band.flux_conversion_factor * 10**(-0.4 * self.B_obs.magnitude)
        flux = self.B_obs.convert_to_flux()
        self.assertEqual(flux.flux, expected_flux)

    def test_convert_to_flux_propagates_Johnson_B_magnitude_uncertainty(self):
        expected_flux = self.B_band.flux_conversion_factor * 10**(-0.4 * self.B_obs.magnitude)
        expected_flux_uncertainty = expected_flux * 0.4 * math.log(10) * self.B_obs.uncertainty
        flux = self.B_obs.convert_to_flux()
        self.assertEqual(flux.flux_uncertainty, expected_flux_uncertainty)

    def test_convert_to_flux_produces_flux_at_band_effective_wavelength(self):
        expected_wavelength = self.B_band.effective_wavelength
        flux = self.B_obs.convert_to_flux()
        self.assertEqual(flux.wavelength, expected_wavelength)

    def test_convert_to_flux_sets_correct_time(self):
        expected_time = 2451663.30
        flux = self.B_obs.convert_to_flux()
        self.assertEqual(flux.time, expected_time)

class TestMagnitudeToFluxConverter(unittest.TestCase):

    def setUp(self):
        self.converter = mag2flux.MagnitudeToFluxConverter()

    def test_calculate_flux_works_for_magnitude_equal_to_zero(self):
        result = self.converter._calculate_flux(0, 1)
        self.assertEqual(1, result)

    def test_calculate_flux_works_for_magnitude_equal_to_one(self):
        result = self.converter._calculate_flux(1, 1)
        self.assertEqual(10**(-0.4), result)

    def test_calculate_flux_works_for_real_data(self):
        expected = 417.5E-11 * 10**(-0.4 * 19.356)
        result = self.converter._calculate_flux(19.356, 417.5E-11)
        self.assertEqual(expected, result)

    def test_calculate_flux_uncertainty_with_zero_mag_uncertainty(self):
        expected = 0
        result = self.converter._calculate_flux_uncertainty(1, 0)
        self.assertEqual(expected, result)

    def test_calculate_flux_uncertainty_with_nonzero_mag_uncertainty(self):
        expected = 0.4 * math.log(10)
        result = self.converter._calculate_flux_uncertainty(1, 1)
        self.assertEqual(expected, result)

    def test_calculate_flux_uncertainty_with_real_data(self):
        flux = 417.5E-11 * 10**(-0.4 * 19.356)
        expected = flux * 0.4 * math.log(10) * 0.02
        result = self.converter._calculate_flux_uncertainty(flux, 0.02)
        self.assertEqual(expected, result)

class TestMagnitudeToFluxConverterConvertMethod(unittest.TestCase):

    def setUp(self):
        self.converter = mag2flux.MagnitudeToFluxConverter()
        self.time = 2451663.30
        self.band = mag2flux.Band('U', 3660.0, 417.5E-11)
        self.observed_mag = mag2flux.ObservedMagnitude(19.356, 
                                                       0.02, 
                                                       self.band, 
                                                       self.time)

    def test_convert_with_real_data_gets_flux_right(self):
        flux = 417.5E-11 * 10**(-0.4 * 19.356)
        expected = mag2flux.MonochromaticFlux(flux, 
                                              flux_uncertainty = 0, 
                                              wavelength = 0, 
                                              time = 0)
        result = self.converter.convert(self.observed_mag)
        self.assertEqual(expected.flux, result.flux)

    def test_convert_with_real_data_gets_flux_uncertainty_right(self):
        flux = 417.5E-11 * 10**(-0.4 * 19.356)
        flux_uncertainty = flux * 0.4 * math.log(10) * 0.02
        expected = mag2flux.MonochromaticFlux(flux, 
                                              flux_uncertainty, 
                                              wavelength = 0, 
                                              time = 0)
        result = self.converter.convert(self.observed_mag)
        self.assertEqual(expected.flux_uncertainty, result.flux_uncertainty)

    def test_convert_with_real_data_does_not_change_wavelength(self):
        expected = mag2flux.MonochromaticFlux(0, 0, self.band.effective_wavelength, self.time)
        result = self.converter.convert(self.observed_mag)
        self.assertEqual(expected.wavelength, result.wavelength)
