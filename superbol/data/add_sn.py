import tables as tb

h5file = tb.open_file('sn_data.h5', 'a')

# This script allows for new supernovae to be added to the SuperBoL
# HDF5 file. Please include your SN data in the fields below as
# indicated in the comments.

sn_name = '' # lowercase SN identifier: i.e. 'sn1987a'

# 1. PHOTOMETRY
# 
# Photometry needs to be added, organized by filter.
# The properties of the filter (effective wavelength, flux zeropoint)
# also need to be added - though default values will be used if none
# are supplied.

filter_name = '' # One-letter filter name such as 'U' or 'g'

# If your know the effective wavelength and flux zeropoint of your
# filter, add them below. If you don't, then leave them set to None
# and default values will be applied based on your filter_name
filter_eff_wl = None
filter_flux_zeropoint = None
