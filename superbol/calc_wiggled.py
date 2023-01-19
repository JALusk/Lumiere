#calc_wiggled.py
import numpy as np
from superbol import flux_wiggler
from superbol import fqbol

num_wiggled_seds = 10

def wiggle_fluxes_n_times(sed):
    wiggled_qbol_fluxes = [None] * num_wiggled_seds
    for i in range(num_wiggled_seds):
        #Make a copy of the original sed
        sed_copy = flux_wiggler.copy_flux_list(sed)
        #Wiggle each flux at each time point in sed
        wiggled_sed = flux_wiggler.wiggle_fluxes(sed_copy)
        #Integrate wiggled sed to get wiggled qbol flux
        wiggled_qbol_fluxes[i] = fqbol.SplineIntegralCalculator().calculate(wiggled_sed)
    return wiggled_qbol_fluxes

def calc_avg_stdev(sed):
    average_qbol_flux = np.average(wiggle_fluxes_n_times(sed))
    print("Average wiggled quasibolometric flux: ", average_qbol_flux)
    stdev_qbol_flux = np.std(wiggle_fluxes_n_times(sed))
    print("STDEV of wiggled quasibolometric fluxes: ", stdev_qbol_flux)
    return [average_qbol_flux, stdev_qbol_flux]
