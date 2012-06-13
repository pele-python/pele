from math import *
import numpy as np
import cpp.ljcpp_ as ljc

class LJ:
    """binary lennard jones potential with smooth cutoff"""
    def __init__(self, eps = 1., sig = 1.):
        print "using lenard jones cpp implementation"
        # self.natoms = natoms
        self.eps = eps
        self.sig = sig

    def getEnergy(self, coords):
        E = ljc.energy(coords, self.eps, self.sig)
        return E

    def getEnergyGradient(self, coords):
        grad=np.zeros(coords.shape[0], np.float64)
        E = ljc.gradient(coords, grad, self.eps, self.sig)
        #print E
        #print coords
        #print grad
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