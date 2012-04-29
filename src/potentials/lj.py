from math import *
import numpy as np #to access np.exp() not built int exp
import potential

def apply_periodic( dr, boxl ):
    for i in xrange(len(dr)):
        #print i, boxl, dr[i],int(dr[i]/(boxl)),
        dr[i] -= boxl*int(dr[i]/(boxl)) 
        #if abs(dr[i]) > (0.5*boxl): dr[i] = boxl - abs(dr[i])
        if dr[i] >= (0.5*boxl): dr[i] -= boxl
        elif dr[i] < (-0.5*boxl): dr[i] += boxl
    #print "dr ",dr
    return dr;

class LJ(potential.potential):
    def __init__(self, eps=1.0, sig=1.0, boxl=None):
        """ simple lennard jones potential"""
        self.sig = sig
        self.eps = eps
        self.boxl = boxl
        if self.boxl == None:
            self.getSep = self.getSep_abs
            self.periodic = False
        else:
            self.getSep = self.getSep_periodic
            self.periodic = True
        print "using Lennard-Jones potential", self.sig, self.eps, 
        if self.periodic: 
            print "periodic with boxl ", self.boxl
        else:
            print ""
        
        self.getEnergy = self.getEnergySlow
        self.getEnergyGradient = self.getEnergyGradientSlow
        if not self.periodic:
            try: 
                import cpp.ljcpp_ as ljc
                print "using fast cpp LJ implimentation"
                self.ljc = ljc
                self.getEnergy = self.getEnergyFast
                self.getEnergyGradient = self.getEnergyGradientFast
            except:
                print "using slow python LJ implimentation"
                pass

    def getSep_periodic(self, vec1, vec2 ):
        assert len(vec1) == 3, "get_sep: vec length not 3"
        dr = (vec2-vec1)
        dr = apply_periodic( dr, self.boxl)
        R = sqrt((dr**2).sum())
        #print R, boxl
        if R > 2./sqrt(3.)*self.boxl:
            print "warning get_sep may not be working right: R = ", R, "boxl=", self.boxl
            print "           ", vec1
            print "           ", vec2
        return R, dr

    def getSep_abs(self, vec1, vec2 ):
        assert len(vec1) == 3, "get_sep: vec length not 3"
        dr = (vec2-vec1)
        R = sqrt((dr**2).sum())
        return R, dr

    def vij(self, r):
        return 4.*self.eps * ( (self.sig/r)**12 - (self.sig/r)**6 )

    def dvij(self, r):
        return 4.*self.eps * ( -12./self.sig*(self.sig/r)**13 + 6./self.sig*(self.sig/r)**7 )

    def getEnergySlow(self, coords):
        natoms = coords.size/3
        V = np.zeros([natoms*3], np.float64) 
        energy=0.
        for i in xrange(natoms):
            for j in xrange(i+1,natoms):
                v1 = coords[i*3:i*3+3]
                v2 = coords[j*3:j*3+3]
                r,dr = self.getSep(v1,v2)
                energy += self.vij(r)
        return energy

    def getEnergyGradientSlow(self, coords):
        natoms = coords.size/3
        V = np.zeros([natoms*3], np.float64) 
        energy=0.
        for i in xrange(natoms):
            for j in xrange(i+1,natoms):
                v1 = coords[i*3:i*3+3]
                v2 = coords[j*3:j*3+3]
                r,dr = self.getSep(v1,v2)
                energy += self.vij(r)
                g = self.dvij(r)
                V[i*3:i*3+3] += -g * dr/r
                V[j*3:j*3+3] += g * dr/r
        return energy,V
    
    def getEnergyFast(self, coords):
        E = self.ljc.energy(coords, self.eps, self.sig)
        return E

    def getEnergyGradientFast(self, coords):
        grad=np.zeros(coords.shape[0], np.float64)
        E = self.ljc.gradient(coords, grad, self.eps, self.sig)
        return E, grad 


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
    

if __name__ == "__main__":
    main()