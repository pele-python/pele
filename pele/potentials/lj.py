from math import *
import numpy as np #to access np.exp() not built int exp

from pygmin.potentials import BasePotential
import fortran.lj as ljf


__all__ = ["LJ"]

class LJ(BasePotential):
    """ simple lennard jones potential"""
    def __init__(self, eps=1.0, sig=1.0, boxl=None):
        self.sig = sig
        self.eps = eps
        self.boxl = boxl
        if self.boxl is None:
            self.periodic = False
            self.boxl = 10000.
        else:
            self.periodic = True

    def getEnergy(self, coords):
        natoms = len(coords) / 3
        E = ljf.ljenergy(
                coords, self.eps, self.sig, self.periodic, self.boxl, [natoms])
        return E

    def getEnergyGradient(self, coords):
        natoms = len(coords) / 3
        E, grad = ljf.ljenergy_gradient(
                coords, self.eps, self.sig, self.periodic, self.boxl, [natoms])
        return E, grad 
    
    def getEnergyList(self, coords, ilist):
        #ilist = ilist_i.getNPilist()
        #ilist += 1 #fortran indexing
        nlist = len(ilist)
        natoms = len(coords) / 3
        E = ljf.energy_ilist(
                coords, self.eps, self.sig, ilist.reshape(-1), self.periodic, 
                self.boxl, [natoms, nlist])
        #ilist -= 1
        return E
    
    def getEnergyGradientList(self, coords, ilist):
        #ilist = ilist_i.getNPilist()
        #ilist += 1 #fortran indexing
        nlist = len(ilist)
        natoms = len(coords) / 3
        E, grad = ljf.energy_gradient_ilist(
                coords, self.eps, self.sig, ilist.reshape(-1), self.periodic, 
                self.boxl, [natoms, nlist])
        #ilist -= 1
        return E, grad 
    
    def getEnergyGradientHessian(self, coords):
        if self.periodic: raise Exception("Hessian not implemented for periodic boundaries")
        from fortran.lj_hess import ljdiff
        g, energy, hess = ljdiff(coords, True, True)
        return energy, g, hess
    



import unittest
class LJTest(unittest.TestCase):
    def setUp(self):
        self.natoms = 10
        self.coords = np.random.uniform(-1,1.,3*self.natoms) * self.natoms**(-1./3)
        self.pot = LJ()
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
    

class TestLJAfterQuench(unittest.TestCase):
    """do the tests after a short quench so that the energies are not crazy large
    """ 
    def setUp(self):
        from pygmin.optimize import mylbfgs
        self.natoms = 10
        self.coords = np.random.uniform(-1,1.,3*self.natoms) * self.natoms**(-1./3)
        self.pot = LJ()
        ret = mylbfgs(self.coords, self.pot, tol=2.)
        self.coords = ret.coords
        self.E = self.pot.getEnergy(self.coords)
        self.Egrad, self.grad = self.pot.getEnergyGradient(self.coords)

    def test_gradient(self):
        e0 = self.pot.getEnergy(self.coords)
        e, g, hess = self.pot.getEnergyGradientHessian(self.coords)
        gnum = self.pot.NumericalDerivative(self.coords)
        maxg = np.max(np.abs(g))
        maxgdiff = np.max(np.abs(g-gnum))
        self.assertAlmostEqual(e0, e, 5)
        self.assertLess(maxgdiff / maxg, 1e-5)

    def test_hessian(self):
        e, g, hess = self.pot.getEnergyGradientHessian(self.coords)
        nhess = self.pot.NumericalHessian(self.coords, eps=1e-8)
#        print "hess", hess
#        print "numerical", nhess
#        diff = hess - nhess
        maxhess = np.max(np.abs(hess))
        maxdiff = np.max(np.abs(hess-nhess))
#        print "diff", hess[:2,:2] - nhess[:2,:2]
#        print maxhess, maxdiff, maxdiff / maxhess
#        print "diff", (hess[:2,:2] - nhess[:2,:2])/maxhess
#        print "diff", (hess - nhess)/maxhess
#        print np.max(np.abs((hess-nhess)/nhess))
        self.assertLess(maxdiff / maxhess, 1e-5)


def main():
    #test class
    natoms = 12
    coords = np.random.uniform(-1,1,natoms*3)*2
    
    lj = LJ()
    E = lj.getEnergy(coords)
    print "E", E 
    E, V = lj.getEnergyGradient(coords)
    print "E", E 
    print "V"
    print V

    print "try a quench"
    from pygmin.optimize import mylbfgs as quench
    quench( coords, lj, iprint=1 )
    #quench( coords, lj.getEnergyGradientNumerical, iprint=1 )

if __name__ == "__main__":
    unittest.main()
    main()
