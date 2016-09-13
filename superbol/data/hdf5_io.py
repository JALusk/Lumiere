import tables as tb

h5file = tb.open_file('sn_data.h5', 'a')

def add_filter_observations(sn_name, filter_name, filter_eff_wl, filter_flux_zeropoint, filter_observations):
    print sn_name
