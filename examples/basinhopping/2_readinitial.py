"""
Example 2: reading coords from file
"""
import numpy as np

from pele.systems import LJCluster

natoms = 12
niter = 100
system = LJCluster(natoms)

coords = np.loadtxt('coords')
coords = coords.reshape(-1)

db = system.create_database()
bh = system.get_basinhopping(database=db)
bh.run(niter)
print "the lowest energy found after", niter, " basinhopping steps is", db.minima()[0].energy
print ""
