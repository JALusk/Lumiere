import unittest

from .context import superbol
from superbol import read_osc
from superbol import lightcurve
from superbol import lqbol

class TestQuasiBolometricLightcurve(unittest.TestCase):

    def setUp(self):
        sn00cb_osc_photometry = read_osc.retrieve_osc_photometry('SN2000cb', path="/home/jlusk/src/superbol/tests/integration/SN2000cb.json")
        fluxes = []
        for photometry_dict in sn00cb_osc_photometry:
            magnitude = read_osc.get_observed_magnitude(photometry_dict)
            fluxes.append(magnitude.convert_to_flux())

        distance = lqbol.Distance(3.0E7 * 3.086E18, 7.0E6 * 3.086E18)
        self.lc_00cb = lightcurve.calculate_lightcurve(fluxes, distance, lqbol.calculate_qbol_luminosity)

    def test_no_negative_luminosities(self):
        self.assertTrue([luminosity.value > 0 for luminosity in self.lc_00cb])
