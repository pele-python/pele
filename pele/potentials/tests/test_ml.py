from __future__ import division
import unittest
import numpy as np
from pele.potentials import MLCost
from pele.optimize import LBFGS_CPP
from numpy.random import randn


def gauss(x, pars):
    mu = pars[0]
    ss = pars[1]
    return np.exp(-(x - mu) ** 2 / (2 * ss)) / np.sqrt(2 * np.pi * ss)


def log_gauss(x, pars):
    mu = pars[0]
    ss = pars[1]
    return (-(x - mu) ** 2 / (2 * ss)) - 0.5 * np.log(2 * np.pi * ss)


class MLTest(unittest.TestCase):
    def setUp(self):
        self.nr_points = 1000
        np.random.seed(42)
        self.mu = 42
        self.ss = 12

    def test_normal_estimates_OK(self):
        observed = np.sqrt(self.ss) * randn(self.nr_points) + self.mu
        pot = MLCost(observed, probf=gauss)
        potl = MLCost(observed, log_probf=log_gauss)
        parameters = [self.mu + 1, self.ss - 2]
        optimizer = LBFGS_CPP(parameters, pot)
        optimizerl = LBFGS_CPP(parameters, potl)
        result = optimizer.run()
        resultl = optimizerl.run()
        opt_parameters = result.coords
        opt_parametersl = resultl.coords
        opt_mu = opt_parameters[0]
        opt_mul = opt_parametersl[0]
        opt_ss = opt_parameters[1]
        opt_ssl = opt_parametersl[1]
        self.assertAlmostEqual(opt_mu, self.mu, delta=2 * self.mu / np.sqrt(self.nr_points))
        self.assertAlmostEqual(opt_ss, self.ss, delta=2 * self.ss / np.sqrt(self.nr_points))
        self.assertAlmostEqual(opt_mu, opt_mul, delta=1e-4)
        self.assertAlmostEqual(opt_ss, opt_ssl, delta=1e-4)
        confidence_intervalsl = potl.get_error_estimate(opt_parametersl)
        for i, par in enumerate(opt_parameters):
            self.assertLessEqual(confidence_intervalsl[i][0], par)
            self.assertLessEqual(par, confidence_intervalsl[i][1])


if __name__ == "__main__":
    unittest.main()

