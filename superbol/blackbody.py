import numpy as np
from scipy.optimize import curve_fit
from astropy import units as u
from astropy import constants as const
from .planck import *


def bb_flux(wavelength, temperature, angular_radius):
    """Observed flux at `wavelength` from a blackbody of `temperature` and `angular_radius` in cgs units

    Args:
        wavelength (float): Wavelength in Angstroms
        temperature (float): Temperature in Kelvin
        angular_radius (float): Angular radius :math:`(\\theta = \\frac{R}{D})`

    Returns:
        Astropy Quantity: blackbody flux in :math:`erg \\; s^{-1} cm^{-2} Anstrom^{-1}`
    """
    bb_flux = (
        (np.pi * u.sr)
        * planck_function(wavelength, temperature)
        * (angular_radius) ** 2
    )

    return bb_flux


def bb_flux_nounits(wavelength, temperature, angular_radius):
    """Identical to bb_flux, but returns the value rather than the Astropy Quantity.

    This function is needed because curve_fit() does not play nice with functions that return Astropy quantities.
    """
    return bb_flux(wavelength, temperature, angular_radius).value






class BlackbodyFit(object):
    def __init__(self):
        self.temperature = None
        self.angular_radius = None
        self.temperature_err = None
        self.angular_radius_err = None
        self.SED = None

    def fit_to_SED(self, SED):
        """Fits the temperature and angular_radius of a blackbody to the SED

        The initial guesses for the `temperature` and `angular_radius` are
        :math:`T = 5000` K and :math:`\\theta = 1.0 \\times 10^{-10}`.
        These are typical for an extragalactic supernovae, but should be used
        with caution for any other objects.

        Args:
            SED (list): List of observed MonochromaticFlux objects.
        """
        wavelengths = [f.wavelength for f in SED]
        fluxes = [f.flux for f in SED]
        flux_uncertainties = [f.flux_uncertainty for f in SED]

        popt, pcov = curve_fit(
            bb_flux_nounits,
            wavelengths,
            fluxes,
            p0=[5000, 1.0e-10],
            sigma=flux_uncertainties,
            absolute_sigma=True,
        )

        self.SED = SED
        self.temperature = popt[0]
        self.angular_radius = popt[1]

        perr = np.sqrt(np.diag(pcov))
        self.temperature_err = perr[0]
        self.angular_radius_err = perr[1]

    def __call__(self, wavelength):
        """Flux of the blackbody at a given wavelength in angstroms

        Args:
            wavelength (float): Wavelength in Angstroms

        Returns:
            bb_flux (float): blackbody flux at given wavelength in cgs units"""
        if self.temperature == None:
            # TODO No test written
            return None
        else:
            bb_flux = (
                np.pi
                * u.sr
                * planck_function(wavelength, self.temperature)
                * self.angular_radius**2
            )
            return bb_flux.value
