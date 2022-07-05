import unittest
import os

from .context import superbol
from astropy.table import Table
from superbol import read_osc
from superbol import lightcurve
from superbol import lqbol
from superbol import lum
from superbol import extinction

dirname = os.path.dirname(__file__)
sn2018hna_extinction = os.path.join(dirname, '../../data/sn2018hna_extinction.dat')
extinction_table = Table.read(sn2018hna_extinction, format = 'ascii')
sn18hna_photometry = os.path.join(dirname, '../../data/SN2018hna.json')

class TestQuasiBolometricLightcurve(unittest.TestCase):

    def setUp(self):
        sn18hna_osc_photometry = read_osc.retrieve_osc_photometry('SN2018hna', path=sn18hna_photometry)
        fluxes = []
        for photometry_dict in sn18hna_osc_photometry:
            try:
                observed_magnitude = read_osc.get_observed_magnitude(photometry_dict)
                extinction_value = extinction.get_extinction_by_name(extinction_table, observed_magnitude.band)
                observed_magnitude.magnitude = extinction.correct_observed_magnitude(observed_magnitude, extinction_value)
                fluxes.append(observed_magnitude.convert_to_flux())
            except:
                pass
        #TODO check distance
        distance = lum.Distance(3.0E7 * 3.086E18, 7.0E6 * 3.086E18)
        self.lc_18hna = lightcurve.calculate_lightcurve(fluxes, distance, lqbol.calculate_qbol_flux)
        #TODO check time values
    def test_18hna_qbol_lightcurve(self):
        print("")
        for luminosity in self.lc_18hna:
            print("{0:9.2f}, {1:5.2f}, 4.1, {2:4.2E}, {3:4.2E}".format(luminosity.time + 2400000.5, luminosity.time + 2400000.5 - 2451656, luminosity.value, luminosity.uncertainty))
