import numpy as np

import NEB.NEB as NEB
import basinhopping
from quench import quench
from potentials.ljcpp import LJ
from mindist.minpermdist_stochastic import minPermDistStochastic as minpermdist
from printing.print_atoms_xyz import printAtomsXYZ


def printpath(fin, coordslist):
    nimages = len(coordslist[:,0])
    for i in range(nimages):
        printAtomsXYZ(fin, coordslist[i,:])


lj = LJ()
natoms = 17
X1 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
ret = quench( X1, lj.getEnergyGradient)
X1 = ret[0]
X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
ret = quench( X2, lj.getEnergyGradient)
X2 = ret[0]

dist, X1, X2 = minpermdist( X1, X2, niter = 10 )
distf = np.linalg.norm(X1 - X2)
print "dist returned        ", dist
print "dist from structures ", distf

#X1 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 0., 1.,] )
#X2 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 1., 0.,] )
import copy
X1i = copy.copy(X1)
X2i = copy.copy(X2)

neb = NEB.NEB(X1, X2, lj )
with open("path.init.xyz", "w") as fout:
    printpath(fout, neb.coords)


neb.optimize()
with open("path.final.xyz", "w") as fout:
    printpath(fout, neb.coords)
