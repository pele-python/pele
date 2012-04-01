from math import *
import numpy as np #to access np.exp() not built int exp

def apply_periodic( dr, boxl ):
    for i in range(len(dr)):
        #print i, boxl, dr[i],int(dr[i]/(boxl)),
        dr[i] -= boxl*int(dr[i]/(boxl)) 
        #if abs(dr[i]) > (0.5*boxl): dr[i] = boxl - abs(dr[i])
        if dr[i] >= (0.5*boxl): dr[i] -= boxl
        elif dr[i] < (-0.5*boxl): dr[i] += boxl
    #print "dr ",dr
    return dr;


class BLJ_interaction_type:
    def __init__(self, eps, sig, rcut):
        self.eps = eps
        self.sig = sig
        self.rcut = rcut*self.sig

        self.ircut2 = 1.0/self.rcut**2
        self.sig6 = self.sig**6
        self.sig12 = self.sig6**2

        sigrc6 = self.sig6/self.rcut**6
        sigrc12 = sigrc6**2

        self.const = 4.0*(sigrc6)-7.0*sigrc12
        self.rconst = (6.0*sigrc12-3.0*sigrc6)/self.rcut**2

class LJpshift:
    """binary lennard jones potential with smooth cutoff"""
    def __init__(self, natoms, ntypeA, boxl=None, rcut=2.5, epsBB=0.5, sigBB=0.88, epsAB=1.5, sigAB=0.8):
        self.boxl = boxl
        self.natoms = natoms
        self.ntypeA = ntypeA


        sigAA = 1.
        epsAA = 1.
        self.AA = BLJ_interaction_type( epsAA, sigAA, rcut)
        self.BB = BLJ_interaction_type( epsBB, sigBB, rcut)
        self.AB = BLJ_interaction_type( epsAB, sigAB, rcut)

        if self.boxl == None:
            self.getSep = self.getSep_abs
            self.periodic = False
        else:
            self.getSep = self.getSep_periodic
            self.periodic = True
        print "using binary Lennard-Jones potential ", self.ntypeA, self.AB.sig, self.BB.sig, self.AB.eps, self.BB.eps
        print "    with cutoff ", rcut,
        if self.periodic: 
            print "periodic with boxl ", self.boxl
        else:
            print ""

    def getSep_periodic(self, vec1, vec2 ):
        assert len(vec1) == 3, "get_sep: vec length not 3"
        dr = (vec2-vec1)
        dr = apply_periodic( dr, self.boxl)
        R2 = ((dr**2).sum())
        return R2, dr

    def getSep_abs( self, vec1, vec2 ):
        assert len(vec1) == 3, "get_sep: vec length not 3"
        dr = (vec2-vec1)
        R2 = ((dr**2).sum())
        return R2, dr

    #define the potential
    #assume R is not more than RCUT
    def vij(self, r2, ir6, T):
        return 4.*T.eps*(T.sig6*ir6*(T.sig6*ir6-1.0) + T.rconst*r2 + T.const)

    #define the first derivative of the potential
    #assume R is not more than RCUT
    #To save some computation, return the gradient / R
    def dvij(self, ir8, ir14, T):
        return -8.0*T.eps*(3.0*(2.0*ir14*(T.sig12)-ir8*T.sig6)-T.rconst)

    def updateEnergy(self, potel, coords, i, j, T):
        r2,dr = self.getSep( coords[range(i*3,i*3+3)], coords[range(j*3,j*3+3)] )
        ir2 = 1./r2
        if ir2 > T.ircut2:
            ir6=ir2**3
            #potel += 4.*T.eps*(T.sig6*ir6*(T.sig6*ir6-1.0) + T.rconst*r2 + T.const)
            potel += self.vij(r2, ir6, T)
        return potel

    def updateEnergyGradient(self, potel, V, coords, i, j, T):
       r2,dr = self.getSep( coords[range(i*3,i*3+3)], coords[range(j*3,j*3+3)] )
       ir2 = 1./r2
       if ir2 > T.ircut2:
           ir6=ir2**3
           #potel += 4.*T.eps*(T.sig6*ir6*(T.sig6*ir6-1.0) + T.rconst*r2 + T.const)
           potel += self.vij(r2, ir6, T)
           ir8=ir6*ir2
           ir14=ir8*ir6
           #g = -8.0*T.eps*(3.0*(2.0*ir14*(T.sig12)-ir8*T.sig6)-T.rconst)
           g = self.dvij(ir8, ir14, T)
           for k in range(3):
              V[i*3+k] += -g * dr[k]
              V[j*3+k] += g * dr[k]
       return potel, V

    def getEnergy(self, coords):
        natoms = self.natoms
        V = np.zeros([natoms*3], np.float64) 
        energy=0.
        for i in range(natoms):
            for j in range(i+1,natoms):
                if i < self.ntypeA and j < self.ntypeA:
                    energy = self.updateEnergy(energy, coords, i,j, self.AA)
                elif i >= self.ntypeA and j >= self.ntypeA:
                    energy = self.updateEnergy(energy, coords, i,j, self.BB)
                else:
                    energy = self.updateEnergy(energy, coords, i,j, self.AB)
        return energy

    def getEnergyGradient(self, coords):
        natoms = self.natoms
        V = np.zeros([natoms*3], np.float64) 
        energy=0.
        for i in range(natoms):
            for j in range(i+1,natoms):
                if i < self.ntypeA and j < self.ntypeA:
                    energy, V = self.updateEnergyGradient(energy, V, coords, i,j, self.AA)
                elif i >= self.ntypeA and j >= self.ntypeA:
                    energy, V = self.updateEnergyGradient(energy, V, coords, i,j, self.BB)
                else:
                    energy, V = self.updateEnergyGradient(energy, V, coords, i,j, self.AB)
                #print "toplot ", r, g, dr[0], dr[1], dr[2], v1, v2
        return energy,V


