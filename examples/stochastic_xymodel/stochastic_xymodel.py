from __future__ import division

import numpy as np
import os

from pele.optimize import StochasticDiagonalLevenbergMarquardt
from pele.potentials import XYModelOnline

class CDFAccumulator(object):
    """
    Compute CDF from data input.
    Convention is that CDF(-\infty) == 1.
    """
    def __init__(self):
        self.total_number = 0
        self.data_x = dict()
    def add_array(self, inp):
        for x in inp:
            self.add(x)
    def add(self, inp):
        already_in = (inp in self.data_x)
        if already_in:
            self.data_x[inp] += 1
        else:
            self.data_x[inp] = 1
        self.total_number += 1
    def get_cdf(self):
        x = self.data_x.keys()
        x = sorted(x)
        cdf_x = []
        remaining_x = self.total_number
        for xi in x:
            cdf_x.append(remaining_x / self.total_number)
            remaining_x -= self.data_x[xi]
        return x, cdf_x
    def get_store_cdf(self):
        self.x, self.cdf_x = self.get_cdf()
        return self.x, self.cdf_x
        
class Ring1D(object):
    def __init__(self, N=8):
        self.N = N
        self.x = np.asarray(range(0, N))
        self.links = []
        self.links.append([0, N - 1])
        for i in xrange(0, N - 1):
            self.links.append([i, i + 1])

if __name__ == "__main__":
    nr_samples = 1000
    nr_spins = 8
    lattice = Ring1D(nr_spins=nr_spins)
    acc = CDFAccumulator()
    pot = XYModelOnline(lattice.nr_spins, lattice.links)
    opt = StochasticDiagonalLevenbergMarquardt(lattice.x, pot)
    for _ in xrange(nr_samples):
        opt.run()
        x_opt = opt.get_x()
        e_opt = pot.get_energy(x_opt)
        acc.add(e_opt)
        opt.reset(lattice.x)
    x, cdf_x = acc.get_cdf()
    print(cdf_x)
