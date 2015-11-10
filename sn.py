import tables as tb
import numpy as np
from astropy import units as u
from mag2flux import mag2flux
from fbol import integrate_fqbol as fqbol_trapezoidal
from fbol import ir_correction, uv_correction_linear, uv_correction_blackbody
from fit_blackbody import bb_fit_parameters
from fit_blackbody import bb_flux_nounits
from specutils import extinction
hdf5_filename = './hdf5/sn_data.h5'

class SN(object):
    """A supernova is the explosion that ends the life of a star

    The SN needs to be conatained within the HDF5 database before it is used
    by SNoBoL. Once there, simply create a supernova by calling the constructor
    with the name of the SN as a string of the form "sn[YEAR][Letter(s)]"

    For example:
    sn1987A = SN('sn1987a')
    sn1999em = SN('sn1999em')

    Attributes
    ----------
    name : Name of the supernova, "sn" followed by the year of first observation
           along with a letter designating the order of observation in that
           year. "sn1987a" was the first SN observed in 1987. "sn2000cb" was the
           eightieth SN observed in 2000.
    """

    def __init__(self, name):
        """Initializes the SN with supplied value for [name]"""
        self.name = name
        self.min_num_obs = 3

        self.read_hdf5()

    def read_hdf5(self):
        h5file = tb.open_file(hdf5_filename, 'r')
        
        self.filter_table = h5file.root.filters
        
        self.sn_node = h5file.get_node('/sn', self.name)
        
        self.get_phot_table()
        self.get_parameter_table()

    def get_phot_table(self):
        self.phot_table = self.sn_node.phot

    def get_parameter_table(self):
        self.parameter_table = self.sn_node.parameters

    def lbol_direct_bh09(self):
        """Calculate the bolometric lightcurve using the direct integration
           method published in Bersten & Hamuy 2009 (2005MNRAS.360..950P)"""
        self.get_observations()
        self.deredden_fluxes()
        self.get_lbol_epochs()
        self.distance_cm, self.distance_cm_err = self.get_distance_cm()
        
        self.lc = np.array([[0.0, 0.0]])
        
        for jd in self.lbol_epochs:
            wavelengths = self.get_wavelengths(jd)
            fluxes = self.get_fluxes(jd)
            flux_errs = self.get_flux_errs(jd)

            fqbol, fqbol_err = fqbol_trapezoidal(wavelengths, fluxes, flux_errs)
            temperature, angular_radius, chisq = bb_fit_parameters(wavelengths,
                                                                   fluxes,
                                                                   flux_errs)
            shortest_wl = np.amin(wavelengths)
            shortest_flux = fluxes[np.argmin(wavelengths)]
            longest_wl = np.amax(wavelengths)

            ir_corr, ir_corr_err = ir_correction(temperature, 
                                                 angular_radius, 
                                                 longest_wl)

#            if shortest_flux < bb_flux_nounits(shortest_wl,
#                                               temperature, 
#                                               angular_radius):
#                uv_corr = uv_correction_linear(shortest_wl, shortest_flux) 
#                uv_corr_err = 0.0 #FIXME
#            else:
#                uv_corr, uv_corr_err = uv_correction_blackbody(temperature,
#                                                               angular_radius,
#                                                               shortest_wl)
            uv_corr, uv_corr_err = uv_correction_blackbody(temperature,
                                                           angular_radius,
                                                           shortest_wl)

            fbol = fqbol + ir_corr + uv_corr
            fbol_err = np.sqrt(np.sum(x*x for x in [fqbol_err, ir_corr_err, uv_corr_err]))
            lum = fbol * 4.0 * np.pi * self.distance_cm**2.0
            self.lc = np.append(self.lc, [[jd, lum]], axis=0)

        self.lc = np.delete(self.lc, (0), axis=0)

    def lbol_bc_bh09_BV(self):
        self.get_magnitudes()
        self.deredden_magnitudes()
        self.get_bc_epochs('B', 'V')
        

    def get_magnitudes(self):
        dtype = [('jd', '>f8'), ('name', 'S1'), ('magnitude', '>f8'), ('uncertainty', '>f8')]
        self.photometry = np.array([(0.0,'0.0',0.0,0.0)], dtype=dtype)
        
        for obs in self.phot_table.iterrows():
            filterid = obs['filter_id']
            for filt in self.filter_table.where('(filter_id == filterid)'):
                self.photometry = np.append(self.photometry, 
                                            np.array([(obs['jd'], 
                                                       filt['name'], 
                                                       obs['magnitude'], 
                                                       obs['uncertainty'])],
                                                     dtype=dtype))

        self.photometry = np.delete(self.photometry, (0), axis=0)

    def deredden_magnitudes(self):
        """Apply the corrections from CCM89 (1989ApJ...345..245C), Table 3"""
        self.Av_gal = self.parameter_table.cols.Av_gal[0]
        self.Av_host = self.parameter_table.cols.Av_host[0]
        self.Av_tot = self.Av_gal + self.Av_host

        ccm89_corr = {'U': 1.569, 'B': 1.337, 'V': 1.0, 'R': 0.751, 'I': 0.479}

        for obs in self.photometry:
            if obs['name'] in ccm89_corr:
                obs['magnitude'] = obs['magnitude'] - ccm89_corr[obs['name']] * self.Av_tot

    def get_distance_cm(self):
        mpc_to_cm = 3.08567758E24
        distance_cm = self.parameter_table.cols.distance_Mpc[0] * mpc_to_cm
        print distance_cm
        distance_cm_err = self.parameter_table.cols.distance_Mpc_err[0] * mpc_to_cm
        return distance_cm, distance_cm_err

    def get_wavelengths(self, jd):
        """Get an array of wavelengths observed on [jd]"""
        wavelengths = np.array([])

        for obs in self.observations:
            if obs[0] == jd:
                wavelengths = np.append(wavelengths, obs[1])

        return wavelengths

    def get_fluxes(self, jd):
        """Get an array of fluxes observed on [jd]"""
        fluxes = np.array([])

        for obs in self.observations:
            if obs[0] == jd:
                fluxes = np.append(fluxes, obs[2])

        return fluxes

    def get_flux_errs(self, jd):
        """Get an array of flux uncertainties observed on [jd]"""
        flux_errs = np.array([])

        for obs in self.observations:
            if obs[0] == jd:
                flux_errs = np.append(flux_errs, obs[3])

        return flux_errs

    def get_lbol_epochs(self):
        """Get only epochs with enough photometric data to calculate Lbol"""
        self.lbol_epochs = np.array([])
        
        for jd_unique in np.unique(self.observations[:,0]):
            num_obs = 0
            for obs in self.observations:
                if obs[0] == jd_unique:
                    num_obs += 1
            if num_obs >= self.min_num_obs:
                self.lbol_epochs = np.append(self.lbol_epochs, jd_unique)

    def get_observations(self):
        self.observations = np.array([[0.0,0.0,0.0,0.0]])
        
        for obs in self.phot_table.iterrows():
            filterid = obs['filter_id']
            for filt in self.filter_table.where('(filter_id == filterid)'):
                flux, flux_err = mag2flux(obs['magnitude'], 
                                          obs['uncertainty'], 
                                          filt['eff_wl'], 
                                          filt['flux_zeropoint'])
                self.observations = np.append(self.observations, [[obs['jd'], 
                                                               filt['eff_wl'], 
                                                               flux, 
                                                               flux_err]],
                                            axis=0)

        self.observations = np.delete(self.observations, (0), axis=0)

    def deredden_fluxes(self):
        self.Av_gal = self.parameter_table.cols.Av_gal[0]
        self.Av_host = self.parameter_table.cols.Av_host[0]
        self.Av_tot = self.Av_gal + self.Av_host

        for obs in self.observations:
            obs[2] = obs[2] * extinction.reddening(obs[1] * u.AA, self.Av_tot, model='ccm89')

    def get_fqbol_trapezoidal(self, jd):
        """Return the quasi-bolometric flux of the supernova on date [jd] using
           trapezoidal integration between the observed fluxes"""

        wavelengths = np.array([])
        fluxes = np.array([])
        flux_err = np.array([])

        for obs in self.observations:
            if obs[0] == jd:
                wavelengths = np.append(wavelengths, obs[1])
                fluxes = np.append(fluxes, obs[2])
                flux_err = np.append(flux_err, obs[3])

        return fqbol_trapezoidal(wavelengths, fluxes, flux_err)
