import numpy as np

import fortran.ljcut as _ljcut
from pygmin.potentials import BasePotential

__all__ = ["LJCut"]

class LJCut(BasePotential):
    """ 
    lennard jones potential with a cutoff that is continuous and smooth
    """
    def __init__(self, eps=1.0, sig=1.0, rcut = 2.5, boxl=None):
        self.sig = sig
        self.eps = eps
        self.rcut = rcut
        self.boxl = boxl
        if self.boxl is None:
            self.periodic = False
            self.boxl = 100000.
        else:
            self.periodic = True
        print "using Lennard-Jones potential", self.sig, self.eps, 
        print "with cutoff", self.rcut,
        if self.periodic: 
            print "periodic with boxl ", self.boxl
        else:
            print ""
        
    def getEnergy(self, coords):
        natoms = len(coords) / 3
        E = _ljcut.ljenergy(
                coords, self.eps, self.sig, self.periodic, self.boxl,
                self.rcut, [natoms])
        return E

    def getEnergyGradient(self, coords):
        natoms = len(coords) / 3
        E, grad = _ljcut.ljenergy_gradient(
                coords, self.eps, self.sig, self.periodic, self.boxl,
                self.rcut, [natoms])
        return E, grad 
    
    def getEnergyList(self, coords, ilist):
        #ilist = ilist_i.getNPilist()
        #ilist += 1 #fortran indexing
        nlist = len(ilist)
        natoms = len(coords) / 3
        E = _ljcut.energy_ilist(
                coords, self.eps, self.sig, ilist.reshape(-1), self.periodic, 
                self.boxl, self.rcut, [natoms, nlist])
        #ilist -= 1
        return E
    
    def getEnergyGradientList(self, coords, ilist):
        #ilist = ilist_i.getNPilist()
        #ilist += 1 #fortran indexing
        nlist = len(ilist)
        natoms = len(coords) / 3
        E, grad = _ljcut.energy_gradient_ilist(
                coords, self.eps, self.sig, ilist.reshape(-1), self.periodic, 
                self.boxl, self.rcut, [natoms, nlist])
        #ilist -= 1
        return E, grad 


import unittest
class LJCutTest(unittest.TestCase):
    def setUp(self):
        self.natoms = 10
        self.coords = np.random.uniform(-1,1.,3*self.natoms) * self.natoms**(-1./3)
        self.pot = LJCut()
        self.E = self.pot.getEnergy(self.coords)
        self.Egrad, self.grad = self.pot.getEnergyGradient(self.coords)
        
        import itertools
        self.ilist = [] #np.zeros([self.natoms*(self.natoms-1)/2, 2])
        k = 0
        for i in range(self.natoms):
            for j in range(i):
                k += 1
                self.ilist.append( [i,j] )
        self.ilist = np.array(self.ilist) 
        #print self.ilist
    
    def test_grad_e(self):
        self.assertAlmostEqual(self.E, self.Egrad, 7)
    def test_lists_e(self):
        e = self.pot.getEnergyList(self.coords, self.ilist)
        self.assertAlmostEqual(self.E, e, 7)
    def test_lists_eg(self):
        e, g = self.pot.getEnergyGradientList(self.coords, self.ilist)
        self.assertAlmostEqual(self.E, e, 7)
        gdiffmax = np.max(np.abs( g-self.grad )) / np.max(np.abs(self.grad))
        self.assertLess(gdiffmax, 1e-7)

if __name__ == "__main__":
    unittest.main()

