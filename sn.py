import tables as tb
import numpy as np
from astropy import units as u
from mag2flux import mag2flux
from fbol import integrate_fqbol as fqbol_trapezoidal
from fbol import ir_correction, uv_correction_linear, uv_correction_blackbody
from fit_blackbody import bb_fit_parameters
from fit_blackbody import bb_flux_nounits
from luminosity import calc_Lbol
from specutils import extinction


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
        self.min_num_obs = 4

        self.read_hdf5()

    def read_hdf5(self):
        """Reads the hdf5 file and returns data on supernova matching [name]"""
        hdf5_filename = './hdf5/sn_data.h5'
        h5file = tb.open_file(hdf5_filename, 'r')
        
        self.filter_table = h5file.root.filters
        
        self.sn_node = h5file.get_node('/sn', self.name)
        
        self.phot_table = self.sn_node.phot
        self.parameter_table = self.sn_node.parameters

    def lbol_direct_bh09(self):
        """Calculate the bolometric lightcurve using the direct integration
           method published in Bersten & Hamuy 2009 (2009ApJ...701..200B)"""
        self.convert_magnitudes_to_fluxes()
        self.deredden_fluxes()
        self.get_lbol_epochs()
        self.distance_cm, self.distance_cm_err = self.get_distance_cm()
        
        self.lc = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
        
        for jd in self.lbol_epochs:
            names = np.array([x['name'] for x in self.converted_obs 
                              if x['jd'] == jd and x['name'] != 'z'])
            wavelengths = np.array([x['wavelength'] for x in self.converted_obs
                                    if x['jd'] == jd and x['name'] != 'z'])
            fluxes = np.array([x['flux'] for x in self.converted_obs
                               if x['jd'] == jd and x['name'] != 'z'])
            flux_errs = np.array([x['uncertainty'] for x in self.converted_obs
                                  if x['jd'] == jd and x['name'] != 'z'])

            sort_indices = np.argsort(wavelengths)
            wavelengths = wavelengths[sort_indices]
            fluxes = fluxes[sort_indices]
            flux_errs = flux_errs[sort_indices]

            fqbol, fqbol_err = fqbol_trapezoidal(wavelengths, fluxes, flux_errs)
            temperature, angular_radius, perr = bb_fit_parameters(wavelengths,
                                                                   fluxes,
                                                                   flux_errs)

            temperature_err = perr[0]
            angular_radius_err = perr[1]

            shortest_wl = np.amin(wavelengths)
            shortest_flux = np.amin(fluxes)
            shortest_flux_err = np.amin(flux_errs)
            longest_wl = np.amax(wavelengths)

            ir_corr, ir_corr_err = ir_correction(temperature,
                                                 temperature_err,
                                                 angular_radius,
                                                 angular_radius_err,
                                                 longest_wl)
            if 'U' in names:
                idx = np.nonzero(names == 'U')[0][0]
                U_flux = fluxes[idx]
                U_wl = wavelengths[idx]
                if U_flux < bb_flux_nounits(U_wl,
                                            temperature, 
                                            angular_radius):
                    uv_corr, uv_corr_err = uv_correction_linear(shortest_wl, 
                                                                shortest_flux, 
                                                                shortest_flux_err) 
                else:
                    uv_corr, uv_corr_err = uv_correction_blackbody(temperature,
                                                             temperature_err,
                                                             angular_radius,
                                                             angular_radius_err,
                                                             shortest_wl)
            else:
                uv_corr, uv_corr_err = uv_correction_blackbody(temperature,
                                                               temperature_err,
                                                               angular_radius,
                                                               angular_radius_err,
                                                               shortest_wl)

            fbol = fqbol + ir_corr + uv_corr
            fbol_err = np.sqrt(np.sum(x*x for x in [fqbol_err, ir_corr_err, uv_corr_err]))
            lum = fbol * 4.0 * np.pi * self.distance_cm**2.0
            lum_err = np.sqrt((4.0 * np.pi * self.distance_cm**2 * fbol_err)**2
                              +(8.0*np.pi * fbol * self.distance_cm * self.distance_cm_err)**2)
            phase = jd - self.parameter_table.cols.explosion_JD[0]
            phase_err = self.parameter_table.cols.explosion_JD_err[0]
            self.lc = np.append(self.lc, [[jd, phase, phase_err, lum, lum_err, temperature, angular_radius]], axis=0)

        self.lc = np.delete(self.lc, (0), axis=0)

        self.write_lbol_plaintext(self.lc, 'direct')

    def lqbol(self):
        """Calculate the quasi-bolometric lightcurve using direct integration
           with trapezoidal integration of the fluxes"""
        self.convert_magnitudes_to_fluxes()
        self.deredden_fluxes()
        self.get_lbol_epochs()
        self.distance_cm, self.distance_cm_err = self.get_distance_cm()
        
        self.qbol_lc = np.array([[0.0, 0.0, 0.0, 0.0, 0.0]])
        
        for jd in self.lbol_epochs:
            wavelengths = np.array([x['wavelength'] for x in self.converted_obs
                                    if x['jd'] == jd and x['name'] != 'z'])
            fluxes = np.array([x['flux'] for x in self.converted_obs
                               if x['jd'] == jd and x['name'] != 'z'])
            flux_errs = np.array([x['uncertainty'] for x in self.converted_obs
                                  if x['jd'] == jd and x['name'] != 'z'])

            sort_indices = np.argsort(wavelengths)
            wavelengths = wavelengths[sort_indices]
            fluxes = fluxes[sort_indices]
            flux_errs = flux_errs[sort_indices]

            fqbol, fqbol_err = fqbol_trapezoidal(wavelengths, fluxes, flux_errs)

            lqbol = fqbol * 4.0 * np.pi * self.distance_cm**2.0
            lqbol_err = np.sqrt((4.0 * np.pi * self.distance_cm**2 * fqbol_err)**2
                              +(8.0*np.pi * fqbol * self.distance_cm * self.distance_cm_err)**2)
            phase = jd - self.parameter_table.cols.explosion_JD[0]
            phase_err = self.parameter_table.cols.explosion_JD_err[0]
            self.qbol_lc = np.append(self.qbol_lc, [[jd, phase, phase_err, lqbol, lqbol_err]], axis=0)
        self.qbol_lc = np.delete(self.qbol_lc, (0), axis=0)
        self.write_lbol_plaintext(self.qbol_lc, 'qbol')

    def lbol_bc_bh09(self, filter1, filter2):
        """Calculate the bolometric lightcurve using the bolometric corrections
           found in Bersten & Hamuy 2009 (2009ApJ...701..200B). These require 
           specifying a color, taken to be filter1 - filter2"""
        self.get_magnitudes()
        self.deredden_UBVRI_magnitudes()
        self.get_bc_epochs(filter1, filter2)
        self.distance_cm, self.distance_cm_err = self.get_distance_cm()

        self.bc_lc = np.array([[0.0, 0.0, 0.0, 0.0, 0.0]])
        
        for i in range(len(self.bc_epochs)):
            jd = self.bc_epochs[i]
            color = self.get_bc_color(jd, filter1, filter2)
            color_err = self.get_bc_color_uncertainty(jd, filter1, filter2) 
            v_mag = np.array([x['magnitude'] for x in self.photometry 
                               if x['jd'] == jd and x['name'] == 'V'])
            v_mag_err = np.array([x['uncertainty'] for x in self.photometry 
                                if x['jd'] == jd and x['name'] == 'V'])      
            lbol_bc, lbol_bc_err = calc_Lbol(color, color_err, filter1+"minus"+filter2, v_mag, v_mag_err, self.distance_cm, self.distance_cm_err)            
            phase = jd - self.parameter_table.cols.explosion_JD[0]
            phase_err = self.parameter_table.cols.explosion_JD_err[0]
            self.bc_lc = np.append(self.bc_lc, [[jd, phase, phase_err, lbol_bc, lbol_bc_err]], axis=0)

        self.bc_lc = np.delete(self.bc_lc, (0), axis=0)
        self.write_lbol_plaintext(self.bc_lc, 'bc_' + filter1 + '-' + filter2)

    def get_bc_color(self, jd, filter1, filter2):
        """Make an array of filter1 - filter 2 on each of the bc_epochs"""

        f1_mag = np.array([x['magnitude'] for x in self.photometry if x['jd'] 
                            == jd and x['name'] == filter1])
        f2_mag = np.array([x['magnitude'] for x in self.photometry if x['jd'] 
                            == jd and x['name'] == filter2])

        return f1_mag - f2_mag

    def get_bc_color_uncertainty(self, jd, filter1, filter2):
        """Make an array of sqrt(dfilter1^2 - dfilter2^2) on each of the 
           bc_epochs"""

        f1_err = np.array([x['uncertainty'] for x in self.photometry if x['jd']
                            == jd and x['name'] == filter1])
        f2_err = np.array([x['uncertainty'] for x in self.photometry if x['jd']
                            == jd and x['name'] == filter2])

        return np.sqrt(f1_err**2 + f2_err**2)

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

    def deredden_UBVRI_magnitudes(self):
        """Apply the corrections from CCM89 (1989ApJ...345..245C), Table 3
        IMPORTANT: This will only deredden the UBVRI magnitudes at the moment"""
        self.Av_gal = self.parameter_table.cols.Av_gal[0]
        self.Av_host = self.parameter_table.cols.Av_host[0]
        self.Av_tot = self.Av_gal + self.Av_host

        ccm89_corr = {'U': 1.569, 'B': 1.337, 'V': 1.0, 'R': 0.751, 'I': 0.479}

        for obs in self.photometry:
            if obs['name'] in ccm89_corr:
                obs['magnitude'] = obs['magnitude'] - ccm89_corr[obs['name']] * self.Av_tot

    def get_bc_epochs(self, filter1, filter2):
        """Get epochs for which observations of both filter1 and filter2 exist"""
        self.bc_epochs = np.array([])
        
        for jd_unique in np.unique(self.photometry['jd']):
            has_filter1 = False
            has_filter2 = False
            for obs in self.photometry:
                if obs['jd'] == jd_unique:
                    if obs['name'] == filter1:
                        has_filter1 = True
                    elif obs['name'] == filter2:
                        has_filter2 = True
            if has_filter1 and has_filter2:
                self.bc_epochs = np.append(self.bc_epochs, jd_unique)

    def get_distance_cm(self):
        mpc_to_cm = 3.08567758E24
        distance_cm = self.parameter_table.cols.distance_Mpc[0] * mpc_to_cm
        distance_cm_err = self.parameter_table.cols.distance_Mpc_err[0] * mpc_to_cm
        return distance_cm, distance_cm_err

    def get_lbol_epochs(self):
        """Get only epochs with enough photometric data to calculate Lbol"""
        self.lbol_epochs = np.array([])
        
        for jd_unique in np.unique(self.converted_obs['jd']):
            num_obs = 0
            for obs in self.converted_obs:
                if obs['jd'] == jd_unique:
                    num_obs += 1
            if num_obs >= self.min_num_obs:
                self.lbol_epochs = np.append(self.lbol_epochs, jd_unique)

    def convert_magnitudes_to_fluxes(self):
        dtype = [('jd', '>f8'), ('name', 'S1'), ('wavelength', '>f8'), ('flux', '>f8'), ('uncertainty', '>f8')]
        self.converted_obs = np.array([(0.0,'0.0',0.0,0.0,0.0)], dtype=dtype)
        
        for obs in self.phot_table.iterrows():
            filterid = obs['filter_id']
            for filt in self.filter_table.where('(filter_id == filterid)'):
                flux, flux_err = mag2flux(obs['magnitude'], 
                                          obs['uncertainty'], 
                                          filt['eff_wl'], 
                                          filt['flux_zeropoint'])
                if 909.09 <= filt['eff_wl'] <= 33333.33:
                    self.converted_obs = np.append(self.converted_obs, 
                                                   np.array([(obs['jd'], 
                                                              filt['name'],
                                                              filt['eff_wl'],
                                                              flux, 
                                                              flux_err)],
                                                            dtype=dtype))

        self.converted_obs = np.delete(self.converted_obs, (0), axis=0)


    def deredden_fluxes(self):
        self.Av_gal = self.parameter_table.cols.Av_gal[0]
        self.Av_host = self.parameter_table.cols.Av_host[0]
        self.Av_tot = self.Av_gal + self.Av_host

        for obs in self.converted_obs:
            obs['flux'] = obs['flux'] * extinction.reddening(obs['wavelength'] * u.AA, self.Av_tot, model='ccm89')

    def write_lbol_plaintext(self, lightcurve, suffix):
        """Write the lightcurve to a file. Append suffix to filename"""

        filename = "lbol_" + self.name + "_" + suffix + ".dat"
        lc_file = open(filename, 'w')
        np.savetxt(lc_file, lightcurve, header = self.name+": JD, Phase, Lbol, err")
        lc_file.close()
