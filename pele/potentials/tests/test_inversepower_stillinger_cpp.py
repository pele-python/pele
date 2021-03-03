from __future__ import division
from __future__ import absolute_import

import unittest
import numpy as np

from pele.potentials import _inversepower_stillinger_cpp
from pele.optimize._quench import lbfgs_cpp

from .test_inversepower_cpp import minimize

class TestInversePowerStillinger_CPP(unittest.TestCase):
    def setUp(self):
        ndim = 2
        npart = 2
        ndof = ndim * npart
        self.x = np.asarray([1.1, 3, 1, 2.1])
        radii = np.asarray([1.05, 1.05])
        exponent = 5
        a = radii[0]*2
        self.etrue = np.power(a / np.sqrt(np.power(self.x[0] - self.x[2], 2) + np.power(self.x[1] - self.x[3], 2)), exponent)
        self.pot = _inversepower_stillinger_cpp.InversePowerStillinger(exponent, radii, ndim=2)
        self.natoms = npart
    def test_energy(self):
        self.assertAlmostEqual(self.etrue, self.pot.getEnergy(self.x))

if __name__ == "__main__":
    unittest.main()

