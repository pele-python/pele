import numpy as np

from pele.potentials import BasePotential
from fortran.soft_sphere_pot import soft_sphere_pot

__all__ = ["SoftSphere"]




class SoftSphere(BasePotential):
    def __init__(self, diams, dimen = 3 ):
        self.diams = diams
        self.dimen = dimen
        
        dmean = np.mean(diams)
        dstd = np.std(diams)
        print "using soft sphere potential with mean diameter", dmean, "st. dev.", dstd
        
    
    def getEnergy(self, coords):
        natoms = len(coords)/self.dimen
        energy, force = soft_sphere_pot(self.dimen, coords, self.diams)
        return energy
    
    def getEnergyGradient(self, coords):
        natoms = len(coords)/self.dimen
        energy, force = soft_sphere_pot(self.dimen, coords, self.diams)
        return energy, force


        
        

def putInBox(coords, boxl):
    natoms = len(coords)/3
    for i in range(len(coords)):
        coords[i] -= boxl*round(coords[i]/boxl)
    return coords
        

def test_soft_sphere(natoms = 9):
    rho = 1.6
    boxl = 1.
    meandiam = boxl / (float(natoms)/rho)**(1./3)
    print "mean diameter", meandiam 
    
    #set up potential
    diams = np.array([meandiam for i in range(natoms)]) #make them all the same
    pot = SoftSphere(diams = diams)

    #initial coordinates
    coords = np.random.uniform(-1,1,[natoms*3]) * (natoms)**(1./3)
    print len(coords)
    E = pot.getEnergy(coords)
    print "initial energy", E 

    printlist = []
    printlist.append((coords.copy(), "intial coords"))
    
    #test a quench with default lbfgs
    from pele.optimize import mylbfgs as quench
    res = quench(coords, pot, iprint=-1)
    printlist.append((res.coords.copy(), "intial coords"))
    print "energy post quench", pot.getEnergy(coords)

    fname = "out.xyz"
    print "saving coordinates to", fname
    from pele.utils.xyz import write_xyz
    with open(fname, "w") as fout:
        for xyz,line2 in printlist:
            xyz = putInBox(xyz, boxl)
            write_xyz(fout, xyz, title=line2) 
            
    
    from scipy.optimize import check_grad
    res = check_grad(pot.getEnergy, pot.getGradient, coords)
    print "testing gradient (should be small)", res

        
        
if __name__ == "__main__":
    test_soft_sphere(120)
