import hdf5_io

# This script allows for new supernovae to be added to the SuperBoL
# HDF5 file. Please include your SN data in the fields below as
# indicated in the comments.

sn_name = 'sn1987a' # lowercase SN identifier: i.e. 'sn1987a'

# 1. PHOTOMETRY
# 
# Photometry needs to be added, organized by filter.
# The properties of the filter (effective wavelength, flux zeropoint)
# also need to be added - though default values will be used if none
# are supplied.

filter_name = '' # One-letter filter name such as 'U' or 'g'

# If your know the effective wavelength and flux zeropoint of your
# filter, add them below. If you don't, then leave them set to None
# and default values will be applied based on your filter_name.
# Flux zeropoints correspond to the flux-at-zero-magnitude in 
# ergs/s/cm2/A
filter_eff_wl = None
filter_flux_zeropoint = None

# Now comes your data. This is a 2D array, each row being:
# JD, magnitude, uncertainty, (note), (reference)

# The full JD of the observation is preferred over MJD
# The note and reference are optional.
# Default values are 'note' and 'ref', respectively.
# If you include a refence for the observation, bibcodes are preferred.

filter_observations = [
        [2450837.8,16.43,0.04,'note','ref'],
        [2450894.7,15.78,0.01,'note','ref'],
        ]

# Once you have entered all of the above information, just uncomment
# the function call and run the script to add those observation to the
# HDF5 file:
hdf5_io.add_filter_observations(sn_name, filter_name, filter_eff_wl, filter_flux_zeropoint, filter_observations)
