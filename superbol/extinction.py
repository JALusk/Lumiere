import numpy as np

def correct_observed_magnitude(observed_magnitude, extinction):
    return observed_magnitude.magnitude - extinction

def get_extinction_by_name(extinction_table, band):
    return extinction_table['A_SandF'][np.where((extinction_table['Filter_name'] == band.name) | (extinction_table['Filter_name'] == band.alt_name))[0][0]]
