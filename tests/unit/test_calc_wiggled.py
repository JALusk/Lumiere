import unittest

from .context import superbol
from superbol import calc_wiggled
from superbol import mag2flux

class TestCalcWiggled(unittest.TestCase):
    
    def setUp(self):
        #runs before each subsequent test in class
        self.flux01 = mag2flux.MonochromaticFlux(100, 2, 1, 0)
        self.flux02 = mag2flux.MonochromaticFlux(200, 2, 2, 0)
        self.flux03 = mag2flux.MonochromaticFlux(150, 2, 3, 0)
        self.sed01 = [self.flux01, self.flux02, self.flux03]

    def test_wiggled_qbol_fluxes_is_list(self):
        result = calc_wiggled.wiggle_fluxes_n_times(self.sed01)
        self.assertEqual(type(result), type([]))

    def test_wiggled_qbol_fluxes_is_n_long(self):
        result = calc_wiggled.wiggle_fluxes_n_times(self.sed01)
        self.assertEqual(len(result), calc_wiggled.num_wiggled_seds)

    def test_average_qbol_flux_above_min(self):
        result = calc_wiggled.calc_avg_stdev(self.sed01)[0]
        min_avg = 339.75 #calculated directly from test values with spline
        self.assertGreaterEqual(result, min_avg)

    def test_average_qbol_flux_below_max(self):
        result = calc_wiggled.calc_avg_stdev(self.sed01)[0]
        max_avg = 347.75 #calculated directly from test values with spline
        self.assertLessEqual(result, max_avg)

    def test_stdev_above_min(self):
        result = calc_wiggled.calc_avg_stdev(self.sed01)[1]
        min_stdev = 0
        self.assertGreaterEqual(result, min_stdev)
    
    def test_stdev_below_max(self):
        result = calc_wiggled.calc_avg_stdev(self.sed01)[1]
        max_stdev = 1600
        self.assertLessEqual(result, max_stdev)
