import numpy as np #to access np.exp() not built int exp
import potential
from scipy import weave
from scipy.weave import converters
from lj import LJ



class ATLJ(potential.potential):
    """
    Lennard Jones + three body Axilrod-Teller term
    
    V = sum_ij VLJ_ij   +  sum_ijk  Z * (1 + 3*cos(t1)*cos(t2*cos(t3)) / (rij * rjk * rik)**3 )
    
    where t1, t2, t3 are the internal angles of the triangle ijk
    
    Z > 0 stabilizes linear vs. triangular geometries 
    """
    def __init__(self, eps=1.0, sig=1.0, Z=1.):
        """ simple lennard jones potential"""
        self.sig = sig
        self.eps = eps
        self.Z = Z
        self.lj = LJ(self.sig, self.eps)

    
    def getEnergy(self, coords):
        """
        use weave inline
        """
        Elj = self.lj.getEnergy(coords)
        
        natoms = coords.size/3
        coords = np.reshape(coords, [natoms,3])
        energy=0.
        Z = self.Z
        #support_code
        code = """
        double drij[3];
        double drik[3];
        double drjk[3];
        energy = 0.;
        for (int i=0; i < natoms; ++i){
            for (int j=0; j<i; ++j){
                for (int k=0; k<j; ++k){
                
                    double rij = 0.;
                    double rik = 0.;
                    double rjk = 0.;
        
                    for (int d=0; d<3; ++d){
                        drij[d] = coords(i,d) - coords(j,d);
                        rij += drij[d]*drij[d];
                    }
                    for (int d=0; d<3; ++d){
                        drjk[d] = coords(j,d) - coords(k,d);
                        rjk += drjk[d]*drjk[d];
                    }
                    for (int d=0; d<3; ++d){
                        drik[d] = coords(i,d) - coords(k,d);
                        rik += drik[d]*drik[d];
                    }
                    
                    rij = sqrt(rij);
                    rjk = sqrt(rjk);
                    rik = sqrt(rik);
                    
                    double ctijk = ( -(drij[0]*drjk[0] + drij[1]*drjk[1] + drij[2]*drjk[2]) / (rij * rjk) );
                    double ctjki = (  (drjk[0]*drik[0] + drjk[1]*drik[1] + drjk[2]*drik[2]) / (rjk * rik) );
                    double ctkij = (  (drik[0]*drij[0] + drik[1]*drij[1] + drik[2]*drij[2]) / (rik * rij) );

                    double r3 = rij*rjk*rik;
                    energy += Z*(1. + 3. * ctijk * ctjki * ctkij) / (r3*r3*r3);
                }
            }
        }
        return_val= energy;
        """
        energy = weave.inline(code, ["coords", "energy", "Z", "natoms"], type_converters=converters.blitz, verbose=2)
        #print "fast energy", Elj, energy
        energy += Elj
        return energy

    def getEnergySlow(self, coords):
        Elj = self.lj.getEnergy(coords)
        
        natoms = coords.size/3
        X = np.reshape(coords, [natoms,3])
        Z = self.Z
        energy = 0.
        for i in range(natoms):
            for j in range(i):
                for k in range(j):
                    #print i, j, k
                    drij = X[i,:] - X[j,:]
                    drik = X[i,:] - X[k,:]
                    drjk = X[j,:] - X[k,:]
                    rij = np.linalg.norm( drij )
                    rik = np.linalg.norm( drik )
                    rjk = np.linalg.norm( drjk )
                    energy += Z * (1. + 3.*\
                            np.dot( drij, -drjk ) * \
                            np.dot(-drij, -drik ) * \
                            np.dot( drik,  drjk ) / (rij*rik*rjk)**2) \
                            / (rij*rik*rjk)**3
        #print "slow energy", Elj, energy
        energy += Elj
        return energy

    

import unittest
class TestATLJ(unittest.TestCase):
    def testenergy(self):
        natoms = 10
        coords = np.random.uniform(-1,1,natoms*3)*2
        
        from optimize.quench import quench
        lj = LJ()
        ret = quench(coords, lj.getEnergyGradient)
        coords = ret[0]
        
        
        atlj = ATLJ(Z=3.)
        e1 = atlj.getEnergy(coords)
        e2 = atlj.getEnergySlow(coords)
        print "%g - %g = %g" % (e1, e2, e1-e2)
        self.assertTrue( abs(e1 - e2) < 1e-12, "ATLJ: two energy methods give different results: %g - %g = %g" % (e1, e2, e1-e2) )



def main():
    #test class
    natoms = 3
    coords = np.random.uniform(-1,1,natoms*3)*2
    
    lj = ATLJ(Z=1.)
    
    
    E = lj.getEnergy(coords)
    print "E", E 
    E, V = lj.getEnergyGradient(coords)
    print "E", E 
    print "V"
    print V

    print "try a quench"
    from optimize.quench import quench
    ret = quench( coords, lj.getEnergyGradient, iprint=-1 )
    #quench( coords, lj.getEnergyGradientNumerical, iprint=1 )
    print "energy ", ret[1]
    print "rms gradient", ret[2]
    print "number of function calls", ret[3]

    from printing.print_atoms_xyz import printAtomsXYZ as printxyz
    coords = ret[0]

    printlist = []
    for i in range(100):
        coords = np.random.uniform(-1,1,natoms*3)*2
        #coords = np.array([0,0,1., 0,0,0, 0,0,2])
        #coords[6:] += np.random.uniform(-1,1,3)*0.1
        ret = quench( coords, lj.getEnergyGradient, iprint=-1 )
        coords = ret[0]
        X = np.reshape(coords, [natoms,3])
        com = X.sum(0) / natoms
        X[:,:] -= com[np.newaxis,:]
        printlist.append(np.reshape(X, natoms*3))
    
    with open("out.xyz", "w") as fout:
        for coords in printlist:
            printxyz(fout, coords)    
    


    
    #print lj.getEnergyGradient(coords)
    #print lj.getEnergyGradientWeave(coords)

if __name__ == "__main__":
    #main()
    unittest.main()

