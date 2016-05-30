# see: test_lbfgs_cpp.py

from __future__ import division

import copy as c
import numpy as np
import os
import unittest as ut

from pele.optimize import StochasticGradientDescent
from pele.optimize import StochasticDiagonalLevenbergMarquardt
from pele.potentials import BasePotential
from pele.potentials import XYModelOnline
from pele.potentials._pele import size_t


class TestSGD_Basic(ut.TestCase):
    def test(self):
        links = np.array([[0, 1], [0, 2], [1, 2]], dtype=size_t)
        x = np.array([np.pi, 0, 0, 0])
        etrue = 2
        pot = XYModelOnline(4, links)
        e = pot.getEnergy(x)
        self.assertAlmostEqual(etrue, e, 10)
        
class TestXY1D8_SGD(ut.TestCase):
    def test(self):
        links = np.array([[0, 7], [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]], dtype=size_t)
        pot = XYModelOnline(8, links)
        # global minimum
        x = np.zeros(8)
        self.assertAlmostEqual(pot.getEnergy(x), -16, 10)
        # local minimum
        x = np.asarray([i * np.pi * 0.25 for i in xrange(8)])
        self.assertAlmostEqual(pot.getEnergy(x), -8 * np.sqrt(2), 10)
        # minimisation
        x = np.asarray(range(8))
        opt1 = StochasticGradientDescent(x, pot, tol=1e-10)
        opt2 = StochasticDiagonalLevenbergMarquardt(x, pot, tol=1e-10)
        opt1.run()
        opt2.run()
        x1 = opt1.get_result().coords
        x2 = opt2.get_result().coords
        self.assertAlmostEqual(pot.getEnergy(x1), -16, 10)
        self.assertAlmostEqual(pot.getEnergy(x2), -16, 10)
        for i in xrange(1, 8):
            d1 = np.abs((x1[i - 1] % (2 * np.pi)) - (x1[i] % (2 * np.pi)))
            d2 = np.abs((x2[i - 1] % (2 * np.pi)) - (x2[i] % (2 * np.pi)))
            self.assertAlmostEqual(d1, 0, 6)
            self.assertAlmostEqual(d1, d2, 6)

if __name__ == "__main__":
    ut.main()

