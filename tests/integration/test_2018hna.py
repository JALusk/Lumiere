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
sn18hna_photometry = os.path.join(dirname, '../../data/SN2018hna.fit')

class TestQuasiBolometricLightcurve(unittest.TestCase):

    def setUp(self):
        sn18hna_fit_photometry = Table.read(sn18hna_photometry, format="fits")
        fluxes = []
        for time, phase, band, magnitude, e_magnitude, telescope, instrument in sn18hna_fit_photometry:
            photometry_dict = {'time' : time, 'band' : band.strip(), 'e_magnitude' : e_magnitude, 'magnitude' : magnitude, 'source': 1}
            try:
                observed_magnitude = read_osc.get_observed_magnitude(photometry_dict)
                extinction_value = extinction.get_extinction_by_name(extinction_table, observed_magnitude.band)
                observed_magnitude.magnitude = extinction.correct_observed_magnitude(observed_magnitude, extinction_value)
                fluxes.append(observed_magnitude.convert_to_flux())
            except:
                pass
        
        distance = lum.Distance(1.282E7 * 3.086E18, 2.02E6 * 3.086E18)
        self.lc_18hna = lightcurve.calculate_lightcurve(fluxes, distance, lqbol.calculate_qbol_flux)

    def test_18hna_qbol_lightcurve(self):
        print("")
        for luminosity in self.lc_18hna:
            print("{0:9.2f}, {1:5.2f}, 4.1, {2:4.2E}, {3:4.2E}".format(luminosity.time, luminosity.time - 2458411.3, luminosity.value, luminosity.uncertainty))




