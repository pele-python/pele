from math import *
import numpy as np #to access np.exp() not built int exp

from pygmin.potentials.potential import potential as basepot
import fortran.lj as ljf


__all__ = ["LJ"]

def apply_periodic( dr, boxl ):
    for i in xrange(len(dr)):
        #print i, boxl, dr[i],int(dr[i]/(boxl)),
        dr[i] -= boxl*int(dr[i]/(boxl)) 
        #if abs(dr[i]) > (0.5*boxl): dr[i] = boxl - abs(dr[i])
        if dr[i] >= (0.5*boxl): dr[i] -= boxl
        elif dr[i] < (-0.5*boxl): dr[i] += boxl
    #print "dr ",dr
    return dr;

class LJ(basepot):
    def __init__(self, eps=1.0, sig=1.0, boxl=None):
        """ simple lennard jones potential"""
        self.sig = sig
        self.eps = eps
        self.boxl = boxl
        if self.boxl is None:
#            self.getSep = self.getSep_abs
            self.periodic = False
            self.boxl = 10000.
        else:
#            self.getSep = self.getSep_periodic
            self.periodic = True
        print "using Lennard-Jones potential", self.sig, self.eps, 
        if self.periodic: 
            print "periodic with boxl ", self.boxl
        else:
            print ""
        
#        self.getEnergy = self.getEnergySlow
#        self.getEnergyGradient = self.getEnergyGradientSlow
#        if not self.periodic:
#            try: 
#                import fortran.lj as ljf
#                print "using fast fortran LJ implementation"
#                self.ljf = ljf
#                self.getEnergy = self.getEnergyFortran
#                self.getEnergyGradient = self.getEnergyGradientFortran
#                self.getEnergyGradientList = self.getEnergyGradientListFortran
#                self.getEnergyList = self.getEnergyListFortran
#            except ImportError:
#                try: 
#                    import cpp.ljcpp_ as ljc
#                    print "using fast cpp LJ implementation"
#                    self.ljc = ljc
#                    self.getEnergy = self.getEnergyFast
#                    self.getEnergyGradient = self.getEnergyGradientFast
#                except ImportError:
#                    print "using slow python LJ implementation"
#            #js850> getEnergyGradientList only implemented for fortran
#            try:
#                import fortran.lj as ljf
#                self.ljf = ljf
#                self.getEnergyGradientList = self.getEnergyGradientListFortran
#            except ImportError:
#                pass
#        else:
#            try: 
#                import fortran.lj as ljf
#                print "using fast fortran LJ implementation"
#                self.ljf = ljf
#                self.getEnergy = self.getEnergyFortran
#                self.getEnergyGradient = self.getEnergyGradientFortran
#                self.getEnergyGradientList = self.getEnergyGradientListFortran
#                self.getEnergyList = self.getEnergyListFortran
#            except ImportError:
#                print "using slow python LJ implementation"

            

#    def getSep_periodic(self, vec1, vec2 ):
#        assert len(vec1) == 3, "get_sep: vec length not 3"
#        dr = (vec2-vec1)
#        dr = apply_periodic( dr, self.boxl)
#        R = sqrt((dr**2).sum())
#        #print R, boxl
#        if R > 2./sqrt(3.)*self.boxl:
#            print "warning get_sep may not be working right: R = ", R, "boxl=", self.boxl
#            print "           ", vec1
#            print "           ", vec2
#        return R, dr
#
#    def getSep_abs(self, vec1, vec2 ):
#        assert len(vec1) == 3, "get_sep: vec length not 3"
#        dr = (vec2-vec1)
#        R = sqrt((dr**2).sum())
#        return R, dr
#
#    def vij(self, r):
#        return 4.*self.eps * ( (self.sig/r)**12 - (self.sig/r)**6 )
#
#    def dvij(self, r):
#        return 4.*self.eps * ( -12./self.sig*(self.sig/r)**13 + 6./self.sig*(self.sig/r)**7 )
#
#    def getEnergySlow(self, coords):
#        natoms = coords.size/3
#        V = np.zeros([natoms*3], np.float64) 
#        energy=0.
#        for i in xrange(natoms):
#            for j in xrange(i+1,natoms):
#                v1 = coords[i*3:i*3+3]
#                v2 = coords[j*3:j*3+3]
#                r,dr = self.getSep(v1,v2)
#                energy += self.vij(r)
#        return energy
#
#    def getEnergyGradientSlow(self, coords):
#        natoms = coords.size/3
#        V = np.zeros([natoms*3], np.float64) 
#        energy=0.
#        for i in xrange(natoms):
#            for j in xrange(i+1,natoms):
#                v1 = coords[i*3:i*3+3]
#                v2 = coords[j*3:j*3+3]
#                r,dr = self.getSep(v1,v2)
#                energy += self.vij(r)
#                g = self.dvij(r)
#                V[i*3:i*3+3] += -g * dr/r
#                V[j*3:j*3+3] += g * dr/r
#        return energy,V
#    
#    def getEnergyFast(self, coords):
#        E = self.ljc.energy(coords, self.eps, self.sig)
#        return E
#
#    def getEnergyGradientFast(self, coords):
#        grad=np.zeros(coords.shape[0], np.float64)
#        E = self.ljc.gradient(coords, grad, self.eps, self.sig)
#        return E, grad 

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
    from pygmin.optimize.quench import quench
    quench( coords, lj.getEnergyGradient, iprint=1 )
    #quench( coords, lj.getEnergyGradientNumerical, iprint=1 )

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
 

if __name__ == "__main__":
    unittest.main()
    main()
