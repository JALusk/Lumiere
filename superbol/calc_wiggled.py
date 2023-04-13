#calc_wiggled.py

import numpy as np
from superbol import flux_wiggler
from superbol import fqbol

import time
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD
from types import FunctionType
from multiprocessing import Pool

num_wiggled_seds = 10

#Parallel via MPI
def wiggle_fluxes_n_times(sed):
    wiggled_qbol_fluxes = [None] * num_wiggled_seds
    pool = Pool() #pulls in MPI infrastructure
    for i in range(num_wiggled_seds):
        #Make a copy of the original sed
        sed_copy = flux_wiggler.copy_flux_list(sed)
        #Wiggle each flux at each time point in sed
        wiggled_sed = flux_wiggler.wiggle_fluxes(sed_copy)
        #Integrate wiggled sed to get wiggled qbol flux
        wiggled_qbol_fluxes[i] = fqbol.SplineIntegralCalculator().calculate(wiggled_sed)
    
    print("Rank ", COMM_WORLD.Get_rank(), "has this list of wiggled fluxes: ", wiggled_qbol_fluxes)
    
    processors = pool.P #P should equal highest rank + 1
    
    return wiggled_qbol_fluxes, processors

#Once at the end - nonparallel
def calc_avg_stdev(sed): 
    pool = Pool()
    #Gather all the wiggled fluxes from each node onto the head node
    all_wiggled_qbol_fluxes = COMM_WORLD.gather(wiggle_fluxes_n_times(sed)[0], root=0)
    #Calculate avg, stdev
    average_qbol_flux = np.average(all_wiggled_qbol_fluxes)
    stdev_qbol_flux = np.std(all_wiggled_qbol_fluxes)

    return [average_qbol_flux, stdev_qbol_flux]

#MPI set-up, see p. 302 in Eff. Comp.
class Pool(object):
    """Process pool using MPI."""
    def __init__(self):
        self.f = None
        self.P = COMM_WORLD.Get_size() #number of processors?
        self.rank = COMM_WORLD.Get_rank()

    def wait(self):
        if self.rank == 0:
            raise RuntimeError("Proc 0 cannot wait!")
        status = MPI.Status()
        while True:
            task = COMM_WORLD.recv(source=0, tag=MPI.ANY_TAG, status=status)
            if not task:
                break
            if isinstance(task, FunctionType):
                self.f = task
                continue
            result = self.f(task)
            COMM_WORLD.isend(result, dest=0, tag=status.tag)

    def map(self, f, tasks):
        N = len(tasks)
        P = self.P #number of worker ranks
        Pless1 = P - 1
        if self.rank != 0:
            self.wait()
            return

        if f is not self.f:
            self.f = f
            requests = []
            for p in range(1, self.P):
                r = COMM_WORLD.isend(f, dest=p)
                requests.append(r)
            MPI.Request.waitall(requests)

        requests = []
        for i, task in enumerate(tasks):
            r = COMM_WORLD.isend(task, dest=(i%Pless1)+1, tag=i)
            requests.append(r)
        MPI.Request.waitall(requests)

        results = []
        for i in range(N):
            result = COMM_WORLD.recv(source=(i%Pless1)+1, tag=i)
            results.append(result)
        return results

    def __del__(self):
        if self.rank == 0:
            for p in range(1, self.P):
                COMM_WORLD.isend(False, dest=p)

#Wiggle fluxes in parallel with MPI and test runtime
def wiggle_in_parallel(sed):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    start = time.time()

    if __name__ == '__main__':
        processors = wiggle_fluxes_n_times(sed)[1] #wiggling in parallel
        print("\nNumber of processors: ", processors)

    calc_avg_stdev(sed)
    
    stop = time.time()
    runtime = stop - start

    print("\nRuntime (sec) for rank ", COMM_WORLD.Get_rank(), " is ", runtime)
    return runtime
