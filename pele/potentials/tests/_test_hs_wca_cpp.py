from __future__ import division
import unittest
import numpy as np
import os
import logging

from pele.optimize import ModifiedFireCPP
from pele.potentials import _hs_wca_cpp
from pele.optimize._quench import lbfgs_cpp
import _base_test

def minimize(coords, pot):
    result = lbfgs_cpp(coords, pot)
    return result.coords, result.energy, result.grad, result.rms

class TestHS_WCA_CPP(_base_test._BaseTest):
    def setUp(self):
        boxlength=1.
        sca = 0.9
        self.natoms = 27 #because the atoms are placed on a plane square box, this integer must have an integer cubic root
        soft_radius= float(boxlength)/(2*np.power(self.natoms,1./3))
        xyz = []
        hs_diameter = float(boxlength)/(np.power(self.natoms,1./3)*(1.+sca))
        hs_radii = np.array([float(hs_diameter/2) for i in xrange(self.natoms)])
        for i in xrange(3):
            for j in xrange(3):
                for k in xrange(3):
                    xyz.append([soft_radius * (1+ 2*i), soft_radius * (1+ 2*j), soft_radius * (1+ 2*k)])
        xyz = np.array(xyz)
        x = xyz.reshape(-1).copy()        
        max_step = 0.99*hs_diameter*sca * .03
        self.xrandom = x + np.random.uniform(-max_step,max_step,[3*self.natoms])
        self.pot = _hs_wca_cpp.HS_WCA(eps=1, sca=sca, radii=hs_radii)
        
        #print "start minimising"
        xyz = minimize(self.xrandom,self.pot)
        #print "end minimising"
        self.xmin = xyz[0].reshape(-1).copy()
        self.Emin = float(xyz[1])
        

class TestErrorPotential(unittest.TestCase):
    def setUp(self):
        self.pot = _hs_wca_cpp._ErrorPotential()
        self.x = np.random.uniform(-1,1,[27*3])
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

# this test was failing regularly, it should be fixed and added back in
#class TestHS_WCA_CPP_NeighborList(_base_test._BaseTest):
#    def setUp(self):
#        self.natoms = 27
#        nlist = [[i,j] for i in xrange(self.natoms) for j in xrange(i+1,self.natoms)]
#        nlist = np.array(nlist, dtype=np.int64).reshape(-1)
#         
#        boxlength=1.
#        sca = 0.9
#        soft_radius= float(boxlength)/(2*np.power(self.natoms,1./3))
#        xyz = []
#        hs_diameter = float(boxlength)/(np.power(self.natoms,1./3)*(1.+sca))
#        hs_radii = np.array([float(hs_diameter/2) for i in xrange(self.natoms)])
#        for i in xrange(3):
#            for j in xrange(3):
#                for k in xrange(3):
#                    xyz.append([soft_radius * (1+ 2*i), soft_radius * (1+ 2*j), soft_radius * (1+ 2*k)])
#        xyz = np.array(xyz)
#        x = xyz.reshape(-1).copy()        
#        max_step = 0.99*hs_diameter*sca * .01
#        self.xrandom = x + np.random.uniform(-max_step,max_step,[3*self.natoms])
#        self.pot = _hs_wca_cpp.HS_WCANeighborList(nlist, eps=1, sca=sca, radii=hs_radii)
#        
#        #print "start minimising"
#        xyz = minimize(self.xrandom,self.pot)
#        #print "end minimising"
#        self.xmin = xyz[0].reshape(-1).copy()
#        self.Emin = float(xyz[1])


if __name__ == "__main__":
    logging.basicConfig(filename='hs_wca_cpp.log',level=logging.DEBUG)
    unittest.main()