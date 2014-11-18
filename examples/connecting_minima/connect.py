"""
example for how to run a double ended connect routine using the system class

We will do the connections for a cluster of 38 Lennard-Jones atoms.
We will load two sets of coordinates from a file, minimize them, and try to find
a connected set of minima and transition states betweeen them.
"""
import numpy as np

from pele.systems import LJCluster

natoms = 38
system = LJCluster(natoms)

# load the coordinates from disk
coords1 = np.genfromtxt("coords.A").flatten()
coords2 = np.genfromtxt("coords.B").flatten()

# minimize them
quencher = system.get_minimizer()
ret1 = quencher(coords1)
ret2 = quencher(coords2)

# make a database and add the minima
db = system.create_database()
coords, E = ret1[:2]
min1 = db.addMinimum(E, coords)
coords, E = ret2[:2]
min2 = db.addMinimum(E, coords)

# do the double ended connect run
connect = system.get_double_ended_connect(min1, min2, db)
connect.connect()
success = connect.success()
if not success:
    print "failed to find connection"
else:
    mints, S, energies = connect.returnPath()
    nts = (len(mints) - 1) / 2
    print "found a path with", nts, "transition states"

    print "plotting energies along the path"
    import pylab as pl

    pl.plot(S, energies, '-')
    pl.xlabel("distance along path")
    pl.ylabel("energy")
    pl.show()
