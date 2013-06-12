import numpy as np
#from potentials.soft_sphere import SoftSphere, putInBox
from pele.potentials.lj import LJ as SoftSphere

np.random.seed(0)

natoms = 120
rho = 1.6
boxl = 1.
meandiam = boxl / (float(natoms)/rho)**(1./3)
print "mean diameter", meandiam 

#set up potential
#diams = np.array([meandiam for i in range(natoms)]) #make them all the same
#pot = SoftSphere(diams = diams)
pot = SoftSphere()


#initial coordinates
coords = np.random.uniform(-1,1,[natoms*3]) * (natoms)**(1./3)
E = pot.getEnergy(coords)
print "initial energy", E 

printlist = [] #list of coordinates saved for printing
printlist.append((coords.copy(), "intial coords"))



#test a quench with default lbfgs
#from optimize.quench import quench
from pele.optimize import lbfgs_ase as quench
coords, E, rms, funcalls = quench(coords, pot.getEnergyGradient, iprint=1)
printlist.append((coords.copy(), "intial coords"))
print "energy post quench", pot.getEnergy(coords)



from scipy.optimize import check_grad
res = check_grad(pot.getEnergy, pot.getGradient, coords)
print "testing gradient (should be small)", res




fname = "out.xyz"
print "saving coordinates to", fname
from pele.printing.print_atoms_xyz import printAtomsXYZ as printxyz
with open(fname, "w") as fout:
    for xyz,line2 in printlist:
        #xyz = putInBox(xyz, boxl)
        printxyz(fout, xyz, line2=line2) 
        
