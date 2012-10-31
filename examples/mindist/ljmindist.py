"""
an example for finding the minimum distance and best alignment
between two lennard jones clusters
"""
import numpy as np

from pygmin.potentials.lj import LJ
from pygmin.optimize import quench
from pygmin.mindist import minPermDistStochastic

pot = LJ()

natoms = 40

#get two random quenched structures to compare
coords1 = np.random.rand(natoms*3)*natoms**(1./3)*1.5
coords2 = np.random.rand(natoms*3)*natoms**(1./3)*1.5
ret1 = quench.lbfgs_py(coords1, pot.getEnergyGradient)
ret2 = quench.lbfgs_py(coords2, pot.getEnergyGradient)
coords1 = ret1[0]
coords2 = ret2[0]

#all the atoms are permutable
permlist = [range(natoms)]

dist, newcoords1, newcoords2 = minPermDistStochastic(coords1, coords2, 
         niter=100, permlist=permlist, verbose=False)

print ""
print "dist =", dist
