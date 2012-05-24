import numpy as np #to access np.exp() not built int exp
import potential
from scipy import weave
from scipy.weave import converters



class LJ(potential.potential):
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

    def getEnergySlow(self, coords):
        natoms = coords.size/3
        coords = np.reshape(coords, [natoms,3])
        energy=0.
        for i in xrange(natoms):
            for j in xrange(i+1,natoms):
                dr = coords[j,:]- coords[i,:]
                r = np.linalg.norm(dr)
                energy += self.vij(r)
        return energy

    def getEnergyGradientSlow(self, coords):
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
    
    def getEnergy(self, coords):
        """
        use weave inline
        """
        natoms = coords.size/3
        energy=0.
        eps = self.eps
        sig2= self.sig**2
        code = """
        double dr[3];
        double r2;
        double g;
        double ir2;
        double ir6;
        double ir12;
        energy = 0.;
        for (int i=0; i < natoms; ++i){
            for (int j=0; j<i; ++j){
                r2 = 0.;
                for (int k=0; k<3; ++k){
                    dr[k] = coords(3*j+k) - coords(3*i+k);
                    r2 += dr[k]*dr[k];
                }
                ir2 = sig2/r2;
                ir6 = ir2*ir2*ir2;
                ir12 = ir6*ir6;
                energy += 4.*eps*( ir12 - ir6 );
                //std::cout << energy << std::endl;
                //g = 4.*eps * ( -12.*ir12 + 6.*ir6 )/r2;
                //for (int k=0; k<3; ++k){
                //    V(3*i+k) += -g * dr[k];
                //    V(3*j+k) +=  g * dr[k];
                //}
            }
        }
        return_val= energy;
        """
        code2 = """
        energy = sig2 + eps;
        """
        energy = weave.inline(code, ["coords", "energy", "sig2", "eps", "natoms"], type_converters=converters.blitz, verbose=2)
        #ret = weave.inline(code2, ["energy", "sig2", "eps"], type_converters=converters.blitz, verbose=2)
        return energy

    
    def getEnergyGradient(self, coords):
        """
        use weave inline
        """
        natoms = coords.size/3
        energy=0.
        V = np.zeros([natoms*3])
        eps = self.eps
        sig2= self.sig**2
        code = """
        double dr[3];
        double r2;
        double g;
        double ir2;
        double ir6;
        double ir12;
        energy = 0.;
        for (int i=0; i < natoms; ++i){
            for (int j=0; j<i; ++j){
                r2 = 0.;
                for (int k=0; k<3; ++k){
                    dr[k] = coords(3*j+k) - coords(3*i+k);
                    r2 += dr[k]*dr[k];
                }
                ir2 = sig2/r2;
                ir6 = ir2*ir2*ir2;
                ir12 = ir6*ir6;
                energy += 4.*eps*( ir12 - ir6 );
                //std::cout << energy << std::endl;
                g = 4.*eps * ( -12.*ir12 + 6.*ir6 )/r2;
                for (int k=0; k<3; ++k){
                    V(3*i+k) += -g * dr[k];
                    V(3*j+k) +=  g * dr[k];
                }
            }
        }
        return_val= energy;
        """
        code2 = """
        energy = sig2 + eps;
        """
        energy = weave.inline(code, ["coords", "V", "energy", "sig2", "eps", "natoms"], type_converters=converters.blitz, verbose=2)
        #ret = weave.inline(code2, ["energy", "sig2", "eps"], type_converters=converters.blitz, verbose=2)
        return energy, V




def main():
    #test class
    natoms = 30
    coords = np.random.uniform(-1,1,natoms*3)*2
    coords0 = coords.copy()
    
    lj = LJ()
    E = lj.getEnergy(coords)
    print "E", E 
    E, V = lj.getEnergyGradientSlow(coords)
    print "E", E 
    print "V"
    print V

    print "try a quench slow"
    from optimize.quench import quench
    ret = quench( coords, lj.getEnergyGradientSlow, iprint=-1 )
    #quench( coords, lj.getEnergyGradientNumerical, iprint=1 )
    print "energy ", ret[1]
    print "rms gradient", ret[2]
    print "number of function calls", ret[3]
    
    print "try a quench weave"
    eweave, vweave = lj.getEnergyGradient(coords)
    ret = quench( coords0, lj.getEnergyGradient, iprint=-1 )
    print "energy weave", eweave
    print "energy weave post quench", ret[1]
    

    
    #print lj.getEnergyGradient(coords)
    #print lj.getEnergyGradientWeave(coords)

if __name__ == "__main__":
    main()
