import unittest
import operator
import os

from .context import superbol
from astropy.table import Table
from superbol import read_osc
from superbol import lightcurve
from superbol import lbc
from superbol import lum
from superbol import extinction

dirname = os.path.dirname(__file__)
sn2000cb_extinction = os.path.join(dirname, '../../data/sn2000cb_extinction.dat')
extinction_table = Table.read(sn2000cb_extinction, format = 'ascii')
sn00cb_photometry = os.path.join(dirname, '../../data/SN2000cb.json')

class TestBolometricCorrectionLightcurve(unittest.TestCase):

    def setUp(self):
        sn00cb_osc_photometry = read_osc.retrieve_osc_photometry('SN2000cb', path=sn00cb_photometry)
        observed_magnitudes = []
        for photometry_dict in sn00cb_osc_photometry:
            try:
                magnitude = read_osc.get_observed_magnitude(photometry_dict)
                extinction_value = extinction.get_extinction_by_name(extinction_table, magnitude.band)
                magnitude.magnitude = extinction.correct_observed_magnitude(magnitude, extinction_value)
                if magnitude.source == "7,13":
                    observed_magnitudes.append(magnitude)
            except:
                pass

        distance = lum.Distance(3.135E7 * 3.086E18, 4.62E6 * 3.086E18)
        self.lc_00cb = lightcurve.calculate_bc_lightcurve(observed_magnitudes, distance, lbc.calculate_bc_flux_h01)

    def test_no_negative_luminosities(self):
        print("")
        for luminosity in self.lc_00cb:
            print("{0:9.2f}, {1:5.2f}, 4.1, {2:4.2E}, {3:4.2E}".format(luminosity.time + 2400000.5, luminosity.time + 2400000.5 - 2451653.8, luminosity.value, luminosity.uncertainty))
