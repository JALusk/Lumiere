import numpy as np
import time
from threading import Thread

#Define structure for threading
class Worker(Thread):
    """Computes x, v, and a of the ith body."""
    def __init__(self, *args, **kwargs):
        super(Worker, self).__init__(*args, **kwargs)
        self.inputs = []
        self.results = []
        self.running = True
        self.daemon = True
        self.start()

    def run(self):
        while self.running:
            if len(self.inputs) == 0:
                continue
            i, x0, v0, G, m, dt = self.inputs.pop(0)
            a_i0 = a(i, x0, G, m)
            v_i1 = a_i0 * dt + v0[i]
            x_i1 = a_i0 * dt**2 + v0[i] * dt + x0[i]
            result = (i, x_i1, v_i1)
            self.results.append(result)

class Pool(object):
    """A collection of P worker threads that distributes tasks
    evenly across them.
    """
    def __init__(self, size):
        self.size = size
        self.workers = [Worker() for p in range(size)]
    def do(self, tasks):
        for p in range(self.size):
            self.workers[p].inputs += tasks[p::self.size]
        while any([len(worker.inputs) != 0 for worker in self.workers]):
            pass
        results = []
        for worker in self.workers:
            results += worker.results
            worker.results.clear()
        return results

    def __del__(self):
        for worker in self.workers:
            worker.running = False

#Define N-body problem
def remove_i(x, i):
    """Drops the ith element of an array."""
    shape = (x.shape[0]-1,) + x.shape[1:]
    y = np.empty(shape, dtype=float)
    y[:i] = x[:i]
    y[i:] = x[i+1:]
    return y

def a(i, x, G, m):
    """The acceleration of the ith mass."""
    x_i = x[i]
    x_j = remove_i(x, i)
    m_j = remove_i(m, i)
    diff = x_j - x_i
    mag3 = np.sum(diff**2, axis=1)**1.5
    result = G * np.sum(diff * (m_j / mag3)[:,np.newaxis], axis=0)
    return result

def timestep(x0, v0, G, m, dt, pool):
    """Computes the next position and velocity for all masses given
    initial conditions and a time step size.
    page 295
    """
    N = len(x0)
    tasks = [(i, x0, v0, G, m, dt) for i in range(N)]
    results = pool.do(tasks)
    x1 = np.empty(x0.shape, dtype=float)
    v1 = np.empty(v0.shape, dtype=float)
    for i, x_i1, v_i1 in results:
        x1[i] = x_i1
        v1[i] = v_i1
    return x1, v1

def initial_cond(N, D):
    """Generates initial conditions for N unity masses at rest
    starting at random positions in D-dimensional space.
    """
    x0 = np.random.rand(N, D)
    v0 = np.zeros((N, D), dtype=float)
    m = np.ones(N, dtype=float)
    return x0, v0, m

#Pool-aware driver function simulate
def simulate(P, N, D, S, G, dt):
    x0, v0, m = initial_cond(N, D)
    pool = Pool(P)
    for s in range(S):
        x1, v1 = timestep(x0, v0, G, m, dt, pool)
        x0, v0 = x1, v1

#Time calculation
Ps = [2, 4, 6, 8]
runtimes = []
i = 0
for P in Ps:
    start = time.time()
    simulate(P, 256, 3, 300, 1.0, 1e-3)
    stop = time.time()
    runtimes.append(stop - start)
    print("Runtime (s) for ", P, "threads is ", runtimes[i])
    i = i + 1

