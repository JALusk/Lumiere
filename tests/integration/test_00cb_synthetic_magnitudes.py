import os
import unittest

from .context import superbol
from astropy.table import Table
from unittest.mock import Mock
from unittest.mock import mock_open
from unittest.mock import patch

from superbol import read_osc
from superbol import lightcurve
from superbol import lqbol
from superbol import lum
from superbol import lbc
from superbol import extinction

dirname = os.path.dirname(__file__)
sn2000cb_extinction = os.path.join(dirname, "../../data/sn2000cb_extinction.dat")
extinction_table = Table.read(sn2000cb_extinction, format="ascii")

sn00cb_synphot_data = """{
"SN2000cb":{
            "photometry":[
                          { 
                           "time":"3.5E41",
                           "band":"U",
                           "e_magnitude":"0.015",
                           "magnitude":"18.906",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"3.5E41",
                           "band":"B",
                           "e_magnitude":"0.015",
                           "magnitude":"17.797",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"3.5E41",
                           "band":"V",
                           "e_magnitude":"0.015",
                           "magnitude":"17.143",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"3.5E41",
                           "band":"R",
                           "e_magnitude":"0.015",
                           "magnitude":"16.778",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"3.5E41",
                           "band":"I",
                           "e_magnitude":"0.015",
                           "magnitude":"16.590",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"3.5E41",
                           "band":"J",
                           "e_magnitude":"0.015",
                           "magnitude":"16.269",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"3.5E41",
                           "band":"H",
                           "e_magnitude":"0.015",
                           "magnitude":"16.166",
                           "u_time":"JD",
		           "source":"1"
                          },
                          { 
                           "time":"4.25E41",
                           "band":"U",
                           "e_magnitude":"0.015",
                           "magnitude":"18.593",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"4.25E41",
                           "band":"B",
                           "e_magnitude":"0.015",
                           "magnitude":"17.564",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"4.25E41",
                           "band":"V",
                           "e_magnitude":"0.015",
                           "magnitude":"16.923",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"4.25E41",
                           "band":"R",
                           "e_magnitude":"0.015",
                           "magnitude":"16.567",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"4.25E41",
                           "band":"I",
                           "e_magnitude":"0.015",
                           "magnitude":"16.366",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"4.25E41",
                           "band":"J",
                           "e_magnitude":"0.015",
                           "magnitude":"16.100",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"4.25E41",
                           "band":"H",
                           "e_magnitude":"0.015",
                           "magnitude":"15.971",
                           "u_time":"JD",
		           "source":"1"
                          },
                          { 
                           "time":"5.00E41",
                           "band":"U",
                           "e_magnitude":"0.015",
                           "magnitude":"18.818",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"5.00E41",
                           "band":"B",
                           "e_magnitude":"0.015",
                           "magnitude":"17.405",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"5.00E41",
                           "band":"V",
                           "e_magnitude":"0.015",
                           "magnitude":"16.679",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"5.00E41",
                           "band":"R",
                           "e_magnitude":"0.015",
                           "magnitude":"16.358",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"5.00E41",
                           "band":"I",
                           "e_magnitude":"0.015",
                           "magnitude":"16.200",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"5.00E41",
                           "band":"J",
                           "e_magnitude":"0.015",
                           "magnitude":"15.872",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"5.00E41",
                           "band":"H",
                           "e_magnitude":"0.015",
                           "magnitude":"15.742",
                           "u_time":"JD",
		           "source":"1"
                          },
                          { 
                           "time":"9.00E41",
                           "band":"U",
                           "e_magnitude":"0.015",
                           "magnitude":"18.225",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"9.00E41",
                           "band":"B",
                           "e_magnitude":"0.015",
                           "magnitude":"17.298",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"9.00E41",
                           "band":"V",
                           "e_magnitude":"0.015",
                           "magnitude":"16.263",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"9.00E41",
                           "band":"R",
                           "e_magnitude":"0.015",
                           "magnitude":"15.696",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"9.00E41",
                           "band":"I",
                           "e_magnitude":"0.015",
                           "magnitude":"14.995",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"9.00E41",
                           "band":"J",
                           "e_magnitude":"0.015",
                           "magnitude":"15.164",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"9.00E41",
                           "band":"H",
                           "e_magnitude":"0.015",
                           "magnitude":"14.998",
                           "u_time":"JD",
		           "source":"1"
                          },
                          { 
                           "time":"10.0E41",
                           "band":"U",
                           "e_magnitude":"0.015",
                           "magnitude":"18.273",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"10.0E41",
                           "band":"B",
                           "e_magnitude":"0.015",
                           "magnitude":"17.321",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"10.0E41",
                           "band":"V",
                           "e_magnitude":"0.015",
                           "magnitude":"16.169",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"10.0E41",
                           "band":"R",
                           "e_magnitude":"0.015",
                           "magnitude":"15.558",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"10.0E41",
                           "band":"I",
                           "e_magnitude":"0.015",
                           "magnitude":"14.953",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"10.00E41",
                           "band":"J",
                           "e_magnitude":"0.015",
                           "magnitude":"15.008",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"10.0E41",
                           "band":"H",
                           "e_magnitude":"0.015",
                           "magnitude":"14.854",
                           "u_time":"JD",
		           "source":"1"
                          },
                          { 
                           "time":"11.0E41",
                           "band":"U",
                           "e_magnitude":"0.015",
                           "magnitude":"17.488",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"11.0E41",
                           "band":"B",
                           "e_magnitude":"0.015",
                           "magnitude":"16.846",
                           "u_time":"JD",
		           "source":"1"
                          },
                          { "time":"11.0E41",
                           "band":"V",
                           "e_magnitude":"0.015",
                           "magnitude":"16.001",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"11.0E41",
                           "band":"R",
                           "e_magnitude":"0.015",
                           "magnitude":"15.497",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"11.0E41",
                           "band":"I",
                           "e_magnitude":"0.015",
                           "magnitude":"14.720",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"11.00E41",
                           "band":"J",
                           "e_magnitude":"0.015",
                           "magnitude":"15.014",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"11.0E41",
                           "band":"H",
                           "e_magnitude":"0.015",
                           "magnitude":"14.871",
                           "u_time":"JD",
		           "source":"1"
                          },
                          { 
                           "time":"12.0E41",
                           "band":"U",
                           "e_magnitude":"0.015",
                           "magnitude":"17.352",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"12.0E41",
                           "band":"B",
                           "e_magnitude":"0.015",
                           "magnitude":"16.800",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"12.0E41",
                           "band":"V",
                           "e_magnitude":"0.015",
                           "magnitude":"15.914",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"12.0E41",
                           "band":"R",
                           "e_magnitude":"0.015",
                           "magnitude":"15.360",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"12.0E41",
                           "band":"I",
                           "e_magnitude":"0.015",
                           "magnitude":"14.627",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"12.00E41",
                           "band":"J",
                           "e_magnitude":"0.015",
                           "magnitude":"14.890",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"12.0E41",
                           "band":"H",
                           "e_magnitude":"0.015",
                           "magnitude":"14.739",
                           "u_time":"JD",
		           "source":"1"
                          },
                          { 
                           "time":"13.0E41",
                           "band":"U",
                           "e_magnitude":"0.015",
                           "magnitude":"16.240",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"13.0E41",
                           "band":"B",
                           "e_magnitude":"0.015",
                           "magnitude":"16.373",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"13.0E41",
                           "band":"V",
                           "e_magnitude":"0.015",
                           "magnitude":"15.821",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"13.0E41",
                           "band":"R",
                           "e_magnitude":"0.015",
                           "magnitude":"15.370",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"13.0E41",
                           "band":"I",
                           "e_magnitude":"0.015",
                           "magnitude":"14.343",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"13.00E41",
                           "band":"J",
                           "e_magnitude":"0.015",
                           "magnitude":"14.922",
                           "u_time":"JD",
		           "source":"1"
                          },
                          {
                           "time":"13.0E41",
                           "band":"H",
                           "e_magnitude":"0.015",
                           "magnitude":"14.788",
                           "u_time":"JD",
		           "source":"1"
                          }
                         ]
           }
}"""


class TestQuasiBolometricLightcurve(unittest.TestCase):
    def setUp(self):
        with patch("builtins.open", mock_open(read_data=sn00cb_synphot_data)):
            self.sn00cb_osc_photometry = read_osc.retrieve_osc_photometry("SN2000cb")
        fluxes = []
        for photometry_dict in self.sn00cb_osc_photometry:
            try:
                observed_magnitude = read_osc.get_observed_magnitude(photometry_dict)
                fluxes.append(observed_magnitude.convert_to_flux())
            except:
                pass

        distance = lum.Distance(3.0e7 * 3.086e18, 7.0e6 * 3.086e18)
        self.lc_00cb = lightcurve.calculate_lightcurve(
            fluxes, distance, lqbol.calculate_qbol_flux
        )

    def test_qbol_lightcurve(self):
        print("")
        for luminosity in self.lc_00cb:
            print(
                "{0:4.2E}, {1:4.2E} +/- {2:4.2E}".format(
                    luminosity.time, luminosity.value, luminosity.uncertainty
                )
            )


class TestBolometricCorrectionLightcurve(unittest.TestCase):
    def setUp(self):
        with patch("builtins.open", mock_open(read_data=sn00cb_synphot_data)):
            sn00cb_osc_photometry = read_osc.retrieve_osc_photometry("SN2000cb")
        observed_magnitudes = []
        for photometry_dict in sn00cb_osc_photometry:
            try:
                magnitude = read_osc.get_observed_magnitude(photometry_dict)
                observed_magnitudes.append(magnitude)
            except:
                pass

        distance = lum.Distance(3.0e7 * 3.086e18, 7.0e6 * 3.086e18)
        self.lc_00cb_bh09 = lightcurve.calculate_bc_lightcurve(
            observed_magnitudes, distance, lbc.calculate_bc_flux_bh09
        )
        self.lc_00cb_h01 = lightcurve.calculate_bc_lightcurve(
            observed_magnitudes, distance, lbc.calculate_bc_flux_h01
        )

    def test_no_negatives_in_bh09_lightcurve(self):
        print("")
        for luminosity in self.lc_00cb_bh09:
            print(
                "{0:4.2E}, {1:4.2E} +/- {2:4.2E}".format(
                    luminosity.time, luminosity.value, luminosity.uncertainty
                )
            )
        self.assertTrue([luminosity.value > 0.0 for luminosity in self.lc_00cb_bh09])

    def test_no_negatives_in_h01_lightcurve(self):
        print("")
        for luminosity in self.lc_00cb_h01:
            print(
                "{0:4.2E}, {1:4.2E} +/- {2:4.2E}".format(
                    luminosity.time, luminosity.value, luminosity.uncertainty
                )
            )
        self.assertTrue([luminosity.value > 0.0 for luminosity in self.lc_00cb_h01])
