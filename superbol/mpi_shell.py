#mpi_shell.py

import numpy as np
from superbol import calc_wiggled
from superbol import mag2flux

import time
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD
from types import FunctionType
from multiprocessing import Pool

#Data from test_calc_wiggled.py
flux01 = mag2flux.MonochromaticFlux(100, 2, 1, 0)
flux02 = mag2flux.MonochromaticFlux(200, 2, 2, 0)
flux03 = mag2flux.MonochromaticFlux(150, 2, 3, 0)
sed = [flux01, flux02, flux03]

#MPI set-up, see p. 302 in Eff. Comp.
class Pool(object):
    """Process pool using MPI."""
    def __init__(self):
        self.f = None
        self.P = COMM_WORLD.Get_size() #number of processors?
        self.rank = COMM_WORLD.Get_rank()

    def getProcessors(self):
        return self.P
    
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



if __name__ == '__main__':
    pool = Pool()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    print("Helloooo")
    #start timer
    start = time.time()
    wiggled_qbol_fluxes = calc_wiggled.wiggle_fluxes_n_times(sed)
    #print("\nNumber of processors: ", processors)
    if rank == 0:
        average, stdev = calc_wiggled.calc_avg_stdev(sed)
        print("Average across all wiggles: ", average)
        print("STDEV across all wiggles: ", stdev)
    stop = time.time()
    runtime = stop - start
    print("\nRuntime (sec) is ", runtime)
