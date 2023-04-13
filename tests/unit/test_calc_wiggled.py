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
        print("result is ", result)
        self.assertEqual(type(result[0]), type([]))

    def test_wiggled_qbol_fluxes_is_n_long(self):
        result = calc_wiggled.wiggle_fluxes_n_times(self.sed01)
        print("result is ", result)
        self.assertEqual(len(result[0]), calc_wiggled.num_wiggled_seds)

    def test_average_qbol_flux_above_min(self):
        result = calc_wiggled.calc_avg_stdev(self.sed01)[0]
        print("Average wiggled quasibolometric flux: ", result)
        min_avg = 339.75 #calculated directly from test values with spline
        self.assertGreaterEqual(result, min_avg)

    def test_average_qbol_flux_below_max(self):
        result = calc_wiggled.calc_avg_stdev(self.sed01)[0]
        max_avg = 347.75 #calculated directly from test values with spline
        self.assertLessEqual(result, max_avg)

    def test_stdev_above_min(self):
        result = calc_wiggled.calc_avg_stdev(self.sed01)[1]
        print("STDEV of wiggled quasibolometric fluxes: ", result)
        min_stdev = 0
        self.assertGreaterEqual(result, min_stdev)
    
    def test_stdev_below_max(self):
        result = calc_wiggled.calc_avg_stdev(self.sed01)[1]
        max_stdev = 1600
        self.assertLessEqual(result, max_stdev)

    def test_wiggling_runtime_above_min(self):
        result = calc_wiggled.wiggle_in_parallel(self.sed01)
        min_runtime = 0
        self.assertGreaterEqual(result, min_runtime)
    
    def test_num_processors_above_min(self):
        result = calc_wiggled.wiggle_fluxes_n_times(self.sed01)[1]
        min_processors = 2
        self.assertGreaterEqual(result, min_processors)
