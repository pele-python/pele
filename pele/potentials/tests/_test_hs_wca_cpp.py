import unittest
import numpy as np
import os
import logging

from pele.potentials import _hs_wca_cpp
from pele.optimize import lbfgs_cpp
import _base_test

def minimize(coords, pot):
    result = lbfgs_cpp(coords, pot)
    return result.coords, result.energy, result.grad, result.rms

class TestHS_WCA_CPP(_base_test._BaseTest):
    def setUp(self):
        sca = 0.99
        a=1./6
        natoms = 9
        xyz = []
        hs_diameters = np.array([1./(3*(1.+sca)) for i in xrange(natoms)])
        for i in xrange(3):
            for j in xrange(3):
                xyz.append([a * (1+ 2*i), a * (1+ 2*j), 0])
        xyz = np.array(xyz)
        x = xyz.reshape(-1).copy()
        self.xrandom = x + np.random.uniform(-1,1,[3*natoms])*0.1
        self.pot = _hs_wca_cpp.HS_WCA(eps=1, sca=sca, radii=(hs_diameters/2), boxl=1)
        
        xyz = minimize(self.xrandom,self.pot)
        self.xmin = xyz[0].reshape(-1).copy()
        self.Emin = float(xyz[1])

class TestErrorPotential(unittest.TestCase):
    def setUp(self):
        self.pot = _hs_wca_cpp._ErrorPotential()
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

class TestHS_WCA_CPP_NeighborList(_base_test._BaseTest):
    def setUp(self):
        self.natoms = 9
        nlist = [[i,j] for i in xrange(self.natoms) for j in xrange(i+1,self.natoms)]
        nlist = np.array(nlist, dtype=np.int64).reshape(-1)
         
        sca = 0.99
        a=1./6
        natoms = self.natoms 
        xyz = []
        hs_diameters = np.array([1./(3*(1.+sca)) for i in xrange(natoms)])
        for i in xrange(3):
            for j in xrange(3):
                xyz.append([a * (1+ 2*i), a * (1+ 2*j), 0])
        xyz = np.array(xyz)
        x = xyz.reshape(-1).copy()
        self.xrandom = x + np.random.uniform(-1,1,[3*natoms])*0.1
        self.pot =_hs_wca_cpp.HS_WCANeighborList(nlist, eps=1, sca=sca, radii=(hs_diameters/2))
        
        xyz = minimize(self.xrandom,self.pot)
        self.xmin = xyz[0].reshape(-1).copy()
        self.Emin = float(xyz[1])



if __name__ == "__main__":
    logging.basicConfig(filename='lj_cpp.log',level=logging.DEBUG)
    unittest.main()