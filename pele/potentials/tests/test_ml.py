from __future__ import division
import unittest
import numpy as np
from pele.potentials import MLCost
from pele.optimize import LBFGS_CPP
from numpy.random import randn

def gauss(x, pars):
    mu = pars[0]
    ss = pars[1]
    return np.exp(-(x - mu)**2 / (2 * ss)) / np.sqrt(2 * np.pi * ss)

class MLTest(unittest.TestCase):
    def setUp(self):
        self.nr_points = 1000
        np.random.seed(42)
        self.mu = 42
        self.ss = 12
    def test_normal_estimates_OK(self):
        observed = np.sqrt(self.ss) * randn(self.nr_points) + self.mu
        pot = MLCost(gauss, observed)
        parameters = [self.mu + 1, self.ss - 2]
        optimizer = LBFGS_CPP(parameters, pot)
        result = optimizer.run()
        opt_parameters = result.coords
        opt_mu = opt_parameters[0]
        print "opt_mu", opt_mu
        opt_ss = opt_parameters[1]
        print "opt_ss", opt_ss
        self.assertAlmostEqual(opt_mu, self.mu, delta=2*self.mu/np.sqrt(self.nr_points))
        self.assertAlmostEqual(opt_ss, self.ss, delta=2*self.ss/np.sqrt(self.nr_points))

if __name__ == "__main__":
    unittest.main()

