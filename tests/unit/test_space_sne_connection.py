import unittest
import requests

from superbol import lightcurve, lqbol, lum, read_osc, read_sne

class TestQuasiBolometricLightcurve(unittest.TestCase):

    def setUp(self):

        self.sn00cb_osc_photometry = read_sne.get_supernova_photometry('SN2000cb')
        fluxes = []
        for photometry_dict in self.sn00cb_osc_photometry:
            try:
                observed_magnitude = read_osc.get_observed_magnitude(
                    photometry_dict)
                fluxes.append(observed_magnitude.convert_to_flux())
            except:
                pass

        distance = lum.Distance(3.0E7 * 3.086E18, 7.0E6 * 3.086E18)
        self.lc_00cb = lightcurve.calculate_lightcurve(
            fluxes, distance, lqbol.calculate_qbol_flux)

    def test_qbol_lightcurve(self):
        print("")
        for luminosity in self.lc_00cb:
            print("{0:4.2E}, {1:4.2E} +/- {2:4.2E}".format(luminosity.time,
                  luminosity.value, luminosity.uncertainty))
