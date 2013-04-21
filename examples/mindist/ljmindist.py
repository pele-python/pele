"""
an example for finding the minimum distance and best alignment
between two lennard jones clusters
"""
import numpy as np

from pygmin.potentials.lj import LJ
from pygmin.optimize._quench_new import lbfgs_py
from pygmin.mindist import MinPermDistAtomicCluster

pot = LJ()

natoms = 40

#get two random quenched structures to compare
coords1 = np.random.rand(natoms*3)*natoms**(1./3)*1.5
coords2 = np.random.rand(natoms*3)*natoms**(1./3)*1.5
ret1 = lbfgs_py(coords1, pot)
ret2 = lbfgs_py(coords2, pot)
coords1 = ret1.coords
coords2 = ret2.coords

#all the atoms are permutable
permlist = [range(natoms)]

mindist = MinPermDistAtomicCluster(niter=100, permlist=permlist, verbose=False)
dist, newcoords1, newcoords2 = mindist(coords1, coords2)

print ""
print "dist =", dist
