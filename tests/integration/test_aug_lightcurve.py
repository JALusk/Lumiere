import unittest

from .context import superbol
from astropy.table import Table
from superbol import read_osc
from superbol import lightcurve
from superbol import lqbol
from superbol import extinction

extinction_table = Table.read("/home/jlusk/src/superbol/data/sn2000cb_extinction.dat", format = 'ascii')

class TestQuasiBolometricLightcurve(unittest.TestCase):

    def setUp(self):
        sn00cb_osc_photometry = read_osc.retrieve_osc_photometry('SN2000cb', path="/home/jlusk/src/superbol/data/SN2000cb.json")
        fluxes = []
        for photometry_dict in sn00cb_osc_photometry:
            try:
                observed_magnitude = read_osc.get_observed_magnitude(photometry_dict)
                extinction_value = extinction.get_extinction_by_name(extinction_table, observed_magnitude.band)
                observed_magnitude.magnitude = extinction.correct_observed_magnitude(observed_magnitude, extinction_value)
                fluxes.append(observed_magnitude.convert_to_flux())
            except:
                pass
        
        distance = lqbol.Distance(3.0E7 * 3.086E18, 7.0E6 * 3.086E18)
        self.lc_00cb = lightcurve.calculate_lightcurve(fluxes, distance, lqbol.calculate_qbol_luminosity)

    def test_no_negative_luminosities(self):
        for luminosity in self.lc_00cb:
            print(luminosity.time, luminosity.time - 51656, 4.1, luminosity.value, luminosity.uncertainty)
        self.assertTrue([luminosity.value > 0 for luminosity in self.lc_00cb])
