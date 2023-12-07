#calc_wiggled.py

import numpy as np
from superbol import flux_wiggler
from superbol import fqbol

import time
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD
from types import FunctionType
# from multiprocessing import Pool
from superbol.mpi_shell import Pool

num_wiggled_seds = 10

#Parallel via MPI
def wiggle_fluxes_n_times(sed):
    wiggled_qbol_fluxes = [None] * num_wiggled_seds
    pool = Pool() #pulls in MPI infrastructure
    if COMM_WORLD.Get_rank() == 0: #Does this only run sometimes?
        for i in range(num_wiggled_seds):
            #Make a copy of the original sed
            sed_copy = flux_wiggler.copy_flux_list(sed)
            #Wiggle each flux at each time point in sed
            wiggled_sed = flux_wiggler.wiggle_fluxes(sed_copy)
            #Integrate wiggled sed to get wiggled qbol flux
            wiggled_qbol_fluxes[i] = fqbol.SplineIntegralCalculator().calculate(wiggled_sed)
    else:
        pool.wait()
    
    num_processors = COMM_WORLD.Get_size()
    return wiggled_qbol_fluxes, num_processors

#Once at the end - nonparallel
def calc_avg_stdev(sed): 
    wiggled_qbol_fluxes = [None] * num_wiggled_seds
    wiggled_qbol_fluxes = wiggle_fluxes_n_times(sed)[0]
    average_qbol_flux = np.average(wiggled_qbol_fluxes)
    print("Average wiggled quasibolometric flux: ", average_qbol_flux)
    stdev_qbol_flux = np.std(wiggled_qbol_fluxes)
    print("STDEV of wiggled quasibolometric fluxes: ", stdev_qbol_flux)
    return [average_qbol_flux, stdev_qbol_flux]

#Wiggle fluxes in parallel with MPI and test runtime
def wiggle_in_parallel(sed):
    start = time.time()

    if __name__ == '__main__':
        processors = wiggle_fluxes_n_times(sed)[1] #wiggling in parallel
        print("\nNumber of processors: ", processors)

    calc_avg_stdev(sed)
    
    stop = time.time()
    runtime = stop - start

    print("\nRuntime (sec) for rank ", COMM_WORLD.Get_rank(), " is ", runtime)
    return runtime