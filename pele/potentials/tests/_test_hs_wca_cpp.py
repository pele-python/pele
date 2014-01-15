from __future__ import division
import unittest
import logging
import numpy as np
import os
import sys

from pele.potentials._hs_wca_cpp import HS_WCA
from pele.utils.xyz import read_xyz
from pele.optimize import lbfgs_cpp
#import _base_test

class TestHS_WCA(unittest.TestCase):
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
        self.xmin = xyz.reshape(-1).copy()
        self.x = self.xmin + np.random.uniform(-1,1,[3*natoms])*0.1
        self.pot =HS_WCA(eps=1, sca=sca, radii=(hs_diameters/2), boxl=1)
        print "Energy= ", self.pot.getEnergy(self.x)
        #result = lbfgs_cpp(xrandom, self.pot, tol=1e-8)
        #self.xmin = result.coords
        #self.xrandom = self.xmin + np.random.uniform(-1,1,[3*natoms])*0.02

    def test_num_grad(self):
        log= logging.getLogger( "BaseTest.grad_test" )
        e, g = self.pot.getEnergyGradient(self.x)
        e1 = self.pot.getEnergy(self.x)
        numerical_g = self.pot.NumericalDerivative(self.x,eps=1e-8)
        log.debug( "g= %r", g )
        log.debug( "numerical_g= %r", numerical_g )
        self.assertLess(np.max(np.abs(g - numerical_g)), 1e-2)
        self.assertAlmostEqual(e, e1, 4)
    
    def test_energy(self):
        e, g = self.pot.getEnergyGradient(self.xmin)    #energy and gradient of minimum configuration generated manually
        e2, g2 = minimize(self.x, self.pot)             #energy and gradient of minimum after minimisation
        self.assertLess(np.max(np.abs(g - g2)), 1e-2)   #assert that gradients are both close to 0
        self.assertAlmostEqual(e, e2, 4)                #assert that energies are both close to 0 
        
def minimize(coords, pot):
    from pele.optimize import lbfgs_cpp
    result = lbfgs_cpp(coords, pot)
    print result.energy, result.rms
    #print result.coords, result.grad
    return result.energy, result.grad

#def main():
#    a=1./6
#    natoms = 9
#    xyz = []
#    diameters = np.array([1./(3*1.4) for i in xrange(natoms)])
#    for i in xrange(3):
#        for j in xrange(3):
#            xyz.append([a * (1+ 2*i), a * (1+ 2*j), 0])
#    xyz = np.array(xyz)
#    xmin = xyz.reshape(-1).copy()
#    xrandom = xmin + np.random.uniform(-1,1,[3*natoms])*0.02
#    pot =HS_WCA(eps=1, sca=0.4, radii=(diameters/2), boxl=1)
#    minimize(xrandom, pot)


if __name__ == "__main__":
    logging.basicConfig(filename='hs_wca.log',level=logging.DEBUG)
    unittest.main()
