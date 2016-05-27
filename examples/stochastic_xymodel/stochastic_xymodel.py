from __future__ import division

import numpy as np
import os

from pele.optimize import LBFGS_CPP
from pele.optimize import SteepestDescentCPP
from pele.optimize import StochasticGradientDescent
from pele.optimize import StochasticDiagonalLevenbergMarquardt
from pele.potentials import XYModelOnline
from pele.potentials._pele import size_t

try:
    import matplotlib.pyplot as p
except Exception as e:
    print(e)

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
    def __init__(self, nr_spins=8, topology="ring"):
        self.nr_spins = nr_spins
        self.x = np.asarray(range(0, nr_spins))
        if topology == "ring":
            self.links = np.zeros((nr_spins, 2), dtype=size_t)
            self.links[0, 0] = 0
            self.links[0, 1] = nr_spins - 1
            for i in xrange(0, nr_spins - 1):
                self.links[i + 1, 0] = i
                self.links[i + 1, 1] = i + 1
        elif topology == "complete":
            self.links = np.zeros((nr_spins * (nr_spins - 1) / 2, 2), dtype=size_t)
            k = 0
            for i in xrange(0, nr_spins):
                for j in xrange(i + 1, nr_spins):
                    self.links[k, :] = np.asarray([i, j])
                    k += 1
        elif topology == "lattice2D":
            L = int(np.sqrt(nr_spins))
            if nr_spins != L * L:
                raise Exception("illegal nr_spins")
            self.links = np.zeros((2 * L * L, 2), dtype=size_t)
            k = 0
            for i in xrange(nr_spins):
                # right
                if (i + 1) % L != 0:
                    # not on right boundary
                    self.links[k, :] = [i, i + 1]
                else:
                    # on right boundary
                    self.links[k, :] = [i, (i + 1) - L]
                k += 1
                # top
                if i // L != (L - 1):
                    # not on top border
                    self.links[k, :] = [i, i + L]
                else:
                    # on top border
                    self.links[k, :] = [i, (i + L) - L * L]
                k += 1
        else:
            raise Exception("illegal topology name: " + topology)

if __name__ == "__main__":
    nr_samples = 100
    nr_spins = 36
    topology = "lattice2D"
    np.random.seed(38)
    lattice = Ring1D(nr_spins=nr_spins, topology=topology)
    acc = CDFAccumulator()
    pot = XYModelOnline(lattice.nr_spins, lattice.links)
    #opt = LBFGS_CPP(lattice.x, pot, tol=1e-15, nsteps=int(1e5))
    #opt = SteepestDescentCPP(lattice.x, pot, tol=1e-10)
    #opt = StochasticGradientDescent(lattice.x, pot, tol=1e-15, nsteps=int(1e5), verbose=False)
    opt = StochasticDiagonalLevenbergMarquardt(lattice.x, pot, tol=1e-15, nsteps=int(1e5))
    for _ in xrange(nr_samples):
        x_ini = 2 * np.pi * np.random.rand(nr_spins)
        opt.reset(x_ini)
        #opt.reset(lattice.x)
        e_ini = pot.getEnergy(x_ini)
        print(e_ini)
        opt.run()
        x_opt = opt.get_result().coords
        e_opt = opt.get_result().energy
        success = opt.get_result().success
        print(success, e_opt, opt.get_result().rms)
        if success:
            acc.add(e_opt)
    x, cdf_x = acc.get_cdf()
    #p.plot(x, cdf_x)
    #p.show()
