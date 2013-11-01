import unittest
import numpy as np
import os

from pele.potentials import _lj_cpp
from pele.utils.xyz import read_xyz
import _base_test

class TestLJ_CPP(_base_test._BaseTest):
    def setUp(self):
        self.pot = _lj_cpp.LJ() 
        self.natoms = 13
        self.xrandom = np.random.uniform(-1,1,[3*self.natoms]) *5.
        current_dir = os.path.dirname(__file__)
        xyz = read_xyz(open(current_dir + "/_lj13_gmin.xyz", "r"))
        self.xmin = xyz.coords.reshape(-1).copy()
        self.Emin = float(xyz.title)

class TestErrorPotential(unittest.TestCase):
    def setUp(self):
        self.pot = _lj_cpp._ErrorPotential()
        self.x = np.random.uniform(-1,1,[9])
    def test(self):
        with self.assertRaises(BaseException):
            self.pot.getEnergy(self.x)
        with self.assertRaises(BaseException):
            self.pot.getEnergyGradient(self.x)
        with self.assertRaises(BaseException):
            self.pot.NumericalDerivative(self.x)
        with self.assertRaises(BaseException):
            self.pot.NumericalHessian(self.x)
#        with self.assertRaises(NotImplementedError):
#            pot.getEnergyGradient(_xrand)

        


if __name__ == "__main__":
    unittest.main()
