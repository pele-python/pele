from __future__ import division

import numpy as np
import os

from pele.optimize import LBFGS_CPP
from pele.optimize import StochasticDiagonalLevenbergMarquardt
from pele.potentials import MLCostOnline

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

from scipy.special import erf

def gauss_cdf(pars, x):
    return 0.5 * (1 - erf((x - pars[0]) / (np.sqrt(2 * pars[1]))))
        
if __name__ == "__main__":
    nr_data_points = 1000
    nr_minimisations = 5
    np.random.seed(38)
    data = np.random.normal(10, 4, nr_data_points)
    x_initial = np.asarray([1., 1.])
    acc = CDFAccumulator()
    pot = MLCostOnline(data, model="gauss")
    #opt = LBFGS_CPP(x_initial, pot, tol=1e-15, nsteps=10000)
    opt = StochasticDiagonalLevenbergMarquardt(x_initial, pot, tol=1e-15, nsteps=10000, verbose=False)
    x_minima = []
    for _ in xrange(nr_minimisations):
        opt.reset(x_initial)
        e_ini = pot.getEnergy(x_initial)
        print(e_ini)
        opt.run()
        x_opt = opt.get_result().coords
        e_opt = opt.get_result().energy
        success = opt.get_result().success
        rms = opt.get_result().rms
        print(success, e_opt, x_opt, rms)
        x_minima.append(x_opt)
    [acc.add(d) for d in data]
    x, cdf_x = acc.get_cdf()
    p.plot(x, cdf_x, label="Data")
    for i, x_opt in enumerate(x_minima):
        p.plot(x, gauss_cdf(x_opt, x), label="Fit " + str(i))
    p.show()
    
