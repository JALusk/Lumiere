import os
import unittest

from .context import superbol
from astropy.table import Table
from superbol import mag2flux
from superbol import extinction

dirname = os.path.dirname(__file__)
sn2000cb_extinction = os.path.join(dirname, '../../data/sn2000cb_extinction.dat')
extinction_table = Table.read(sn2000cb_extinction, format = 'ascii')

class TestExtinctionCorrection(unittest.TestCase):

    def setUp(self):
        self.B_band = mag2flux.Band('B', 'CTIO B', 4380.0, 632.0E-11)
        self.V_band = mag2flux.Band('V', 'CTIO V', 5450.0, 363.1E-11)
        self.B_obs = mag2flux.ObservedMagnitude(18.793, 0.02, self.B_band,
                                                2451663.30)
        self.V_obs = mag2flux.ObservedMagnitude(18.078, 0.015, self.V_band,
                                                2451663.30)
        self.A_B = 0.409
        self.A_V = 0.302
        self.sn_name = 'SN2000cb'

    def test_correct_observed_magnitude(self):
        expected = self.B_obs.magnitude - self.A_B
        result  = extinction.correct_observed_magnitude(self.B_obs, self.A_B)
        self.assertEqual(expected, result)

    def test_get_B_extinction_by_name(self):
        expected = 0.409
        result = extinction.get_extinction_by_name(extinction_table, self.B_band)
        self.assertEqual(expected, result)

    def test_get_V_extinction_by_name(self):
        expected = 0.302
        result = extinction.get_extinction_by_name(extinction_table, self.V_band)
        self.assertEqual(expected, result)


