import unittest
import numpy as np
import os
import logging

from pele.potentials import _wca_cpp
from pele.optimize._quench import lbfgs_cpp, modifiedfire_cpp
import _base_test

def minimize(coords, pot):
    result = lbfgs_cpp(coords, pot)
    return result.coords, result.energy, result.grad, result.rms

class TestWCA_CPP(_base_test._BaseTest):
    def setUp(self):
        boxv=[3,3,3]
        self.pot = _wca_cpp.WCA(ndim=3, boxvec = boxv) 
        self.natoms = 25
        self.xrandom = np.random.uniform(-1,1,[3*self.natoms])*0.001
        
        xyz = np.random.uniform(-1,1,[3*self.natoms])
        xyz = minimize(xyz,self.pot)
        xyz = minimize(xyz[0],self.pot)
        self.xmin = xyz[0].reshape(-1).copy()
        self.Emin = float(xyz[1])

class TestWCA2D_CPP(_base_test._BaseTest):
    def setUp(self):
        boxv=[3,3]
        self.pot = _wca_cpp.WCA(ndim=2,boxvec = boxv) 
        self.natoms = 25
        self.xrandom = np.random.uniform(-1,1,[2*self.natoms])*0.01
        xy = np.random.uniform(-1,1,[2*self.natoms])
        xy = minimize(xy,self.pot)
        xy = minimize(xy[0],self.pot)
        self.xmin = xy[0].reshape(-1).copy()
        self.Emin = float(xy[1])

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
