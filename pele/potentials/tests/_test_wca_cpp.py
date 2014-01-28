import unittest
import numpy as np
import os
import logging

from pele.potentials import _wca_cpp
from pele.optimize import lbfgs_cpp
import _base_test

def minimize(coords, pot):
    result = lbfgs_cpp(coords, pot)
    return result.coords, result.energy, result.grad, result.rms

class TestWCA_CPP(_base_test._BaseTest):
    def setUp(self):
        self.pot = _wca_cpp.WCA() 
        self.natoms = 13
        self.xrandom = np.random.uniform(-1,1,[3*self.natoms]) *5.
        xyz = minimize(self.xrandom,self.pot)
        self.xmin = xyz[0].reshape(-1).copy()
        self.Emin = float(xyz[1])

class TestErrorPotential(unittest.TestCase):
    def setUp(self):
        self.pot = _wca_cpp._ErrorPotential()
        self.x = np.random.uniform(-1,1,[9])
    def test(self):
        with self.assertRaises(RuntimeError):
            self.pot.getEnergy(self.x)
        with self.assertRaises(RuntimeError):
            self.pot.getEnergyGradient(self.x)
        with self.assertRaises(RuntimeError):
            self.pot.NumericalDerivative(self.x)
        with self.assertRaises(RuntimeError):
            self.pot.NumericalHessian(self.x)
#        with self.assertRaises(NotImplementedError):
#            pot.getEnergyGradient(_xrand)

class TestWCA_CPP_NeighborList(_base_test._BaseTest):
    def setUp(self):
        self.natoms = 13
        nlist = [[i,j] for i in xrange(self.natoms) for j in xrange(i+1,self.natoms)]
        nlist = np.array(nlist, dtype=np.int64).reshape(-1)
        self.pot = _wca_cpp.WCANeighborList(nlist) 
        self.xrandom = np.random.uniform(-1,1,[3*self.natoms]) *5.
        xyz = minimize(self.xrandom,self.pot)
        self.xmin = xyz[0].reshape(-1).copy()
        self.Emin = float(xyz[1])



if __name__ == "__main__":
    logging.basicConfig(filename='lj_cpp.log',level=logging.DEBUG)
    unittest.main()
