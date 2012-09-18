import numpy as np

import fortran.ljcut as _ljcut
from pygmin.potentials.potential import potential as basepot

class LJCut(basepot):
    def __init__(self, eps=1.0, sig=1.0, rcut = 2.5, boxl=None):
        """ 
        simple lennard jones potential with a cutoff that is continuous and smooth
        """
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
        ilist += 1 #fortran indexing
        nlist = len(ilist)
        natoms = len(coords) / 3
        E = _ljcut.energy_ilist(
                coords, self.eps, self.sig, ilist, self.periodic, 
                self.boxl, self.rcut, [natoms, nlist])
        ilist -= 1
        return E
    
    def getEnergyGradientList(self, coords, ilist):
        #ilist = ilist_i.getNPilist()
        ilist += 1 #fortran indexing
        nlist = len(ilist)
        natoms = len(coords) / 3
        E, grad = _ljcut.energy_gradient_ilist(
                coords, self.eps, self.sig, ilist, self.periodic, 
                self.boxl, self.rcut, [natoms, nlist])
        ilist -= 1
        return E, grad 
