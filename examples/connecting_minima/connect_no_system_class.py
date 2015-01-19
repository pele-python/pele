"""
example for how to run a double ended connect routine without using the system class.
We strongly recommend not doing it this way and setting up a system class first.  It
will be much easier for you.

We will do the connections for a cluster of 38 Lennard-Jones atoms.
We will load two sets of coordinates from a file, minimize them, and try to find
a connected set of minima and transition states betweeen them.
"""

import numpy as np

from pele.potentials.lj import LJ
from pele.optimize import lbfgs_py
from pele.landscape import DoubleEndedConnect, smooth_path
from pele.mindist import MinPermDistAtomicCluster
from pele.transition_states import orthogopt
from pele.storage import Database
from pele.utils.xyz import write_xyz

np.random.seed(0)

# set up the potential
pot = LJ()

# import the starting and ending points and quench them, 
coords1 = np.genfromtxt("coords.A")
coords2 = np.genfromtxt("coords.B")
res1 = lbfgs_py(coords1.reshape(-1), pot)
res2 = lbfgs_py(coords2.reshape(-1), pot)
coords1 = res1.coords
coords2 = res2.coords
E1 = res1.energy
E2 = res2.energy
natoms = len(coords1) / 3

# add the minima to a database
dbfile = "database.sqlite"
database = Database(dbfile)
database.addMinimum(E1, coords1)
database.addMinimum(E2, coords2)
min1 = database.minima()[0]
min2 = database.minima()[1]




# set up the structural alignment routine.
# we have to deal with global translational, global rotational,
# and permutational symmetry.
permlist = [range(natoms)]
mindist = MinPermDistAtomicCluster(permlist=permlist, niter=10)

# The transition state search needs to know what the eigenvector corresponding
# to the lowest nonzero eigenvector is.  For this we need to know what the
# eivenvector corresponding to the zero eigenvalues are.  These are related
# to global symmetries.  for this system we have 3 zero eigenvalues for translational
# symmetry and 3 for rotational symmetry.  The function orthogopt
# will take care of this for us.
orthogZeroEigs = orthogopt
tsSearchParams = dict()
tsSearchParams["orthogZeroEigs"] = orthogZeroEigs
tsSearchParams["iprint"] = 10

if True:
    # set additional parameters to make it run faster
    tsSearchParams["nsteps_tangent1"] = 3
    tsSearchParams["nsteps_tangent2"] = 25
    tsSearchParams["lowestEigenvectorQuenchParams"] = {"nsteps": 50}
    tsSearchParams["nfail_max"] = 200

    NEBquenchParams = dict()
    NEBquenchParams["iprint"] = 100

    NEBparams = dict()
    NEBparams["k"] = 5.
    NEBparams["image_density"] = 15
    NEBparams["NEBquenchParams"] = NEBquenchParams

    local_connect_params = dict()
    local_connect_params["tsSearchParams"] = tsSearchParams
    local_connect_params["NEBparams"] = NEBparams

myconnect = DoubleEndedConnect(min1, min2, pot, mindist, database,
                               local_connect_params=local_connect_params,
)

if (myconnect.graph.areConnected(min1, min2)):
    print "ALERT: The minima are already connected in the database file", dbfile
    print "       Delete it for a fresh run."

myconnect.connect()
print ""
print "found a path!"



# the path is now saved in the database.  Lets retrieve it and 
# print it in a more visual format
# now retrieve the path for printing
print ""
mints, S, energies = myconnect.returnPath()
nmin = (len(mints) - 1) / 2 + 1
nts = nmin - 1
print "the path has %d minima and %d transition states" % (nmin, nts)
eofs = "path.EofS"
print "saving energies to", eofs
with open(eofs, "w") as fout:
    for i in range(len(S)):
        fout.write("%f %f\n" % (S[i], energies[i]))

xyzfile = "path.xyz"
print "saving path in xyz format to", xyzfile
with open(xyzfile, "w") as fout:
    for m in mints:
        write_xyz(fout, m.coords, title=str(m.energy))

xyzfile = "path.smooth.xyz"
print "saving smoothed path in xyz format to", xyzfile
clist = [m.coords for m in mints]
smoothed = smooth_path(clist, mindist)
with open(xyzfile, "w") as fout:
    for coords in smoothed:
        write_xyz(fout, coords)

if False:
    try:
        import matplotlib.pyplot as plt

        plt.plot(S, energies)
        plt.show()
    except ImportError:
        print "problem plotting with pyplot, skipping"
