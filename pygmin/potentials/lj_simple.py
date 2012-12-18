import numpy as np #to access np.exp() not built int exp

from pygmin.potentials import BasePotential


class LJ(BasePotential):
    """
    A stripped down, simple version of LJ for example purposes
    """
    def __init__(self, eps=1.0, sig=1.0):
        """ simple lennard jones potential"""
        self.sig = sig
        self.eps = eps

    def vij(self, r):
        return 4.*self.eps * ( (self.sig/r)**12 - (self.sig/r)**6 )

    def dvij(self, r):
        return 4.*self.eps * ( -12./self.sig*(self.sig/r)**13 + 6./self.sig*(self.sig/r)**7 )

    def getEnergy(self, coords):
        natoms = coords.size/3
        coords = np.reshape(coords, [natoms,3])
        energy=0.
        for i in xrange(natoms):
            for j in xrange(i+1,natoms):
                dr = coords[j,:]- coords[i,:]
                r = np.linalg.norm(dr)
                energy += self.vij(r)
        return energy

    def getEnergyGradient(self, coords):
        natoms = coords.size/3
        coords = np.reshape(coords, [natoms,3])
        energy=0.
        V = np.zeros([natoms,3])
        for i in xrange(natoms):
            for j in xrange(i+1,natoms):
                dr = coords[j,:]- coords[i,:]
                r = np.linalg.norm(dr)
                energy += self.vij(r)
                g = self.dvij(r)
                V[i,:] += -g * dr/r
                V[j,:] += g * dr/r
        V = V.reshape([natoms*3])
        return energy,V



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
    ret = quench( coords, lj.getEnergyGradient, iprint=-1 )
    #quench( coords, lj.getEnergyGradientNumerical, iprint=1 )
    print "energy ", ret[1]
    print "rms gradient", ret[2]
    print "number of function calls", ret[3]
        

if __name__ == "__main__":
    main()
