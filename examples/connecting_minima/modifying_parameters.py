"""
example for how to modify parameters in the double ended connect routine

This example will be based on a cluster of 38 Lennard-Jones atoms.
We will load two sets of coordinates from a file, minimize them, and try to find
a connected set of minima and transition states betweeen them.

For more information about how to change parameters, see the documentation
for the system class, the parameter tree and DoubleEndedConnect.

See the class LJCluster for what the default parameters are for this system
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


#
# here is where we modify the parameters
#
"""
double ended connect is a heierarchical routine where you have functions calling functions which
then call functions.  In order to reproduce this logically the parameters are stored in 
heierarchichal dictionaries in the system class starting with system.params 

there are many more parameters which can be adjusted.  Some can have a large effect
on how fast you find a connection and how direct that connection is.  There are
additional tools in the gui (neb explorer, connect explorer, etc) for finding good 
parameters.

Ultimately, if you find a good parameter set for your system, you should add the
changes to your system class.  Then they will be set by default every time you
call system.get_double_ended_connect()
"""
# specify a fresh_connect.  i.e. ignore all existing minima in the database
decparams = system.params.double_ended_connect
decparams.fresh_connect = True
# specify the maximum number of NEB transition state candidates to refine
decparams.local_connect_params.nrefine_max = 5

# NEB parameters:
# reinterpolate the path to achieve equidistant spacing every 20 steps
NEBparams = decparams.local_connect_params.NEBparams
NEBparams.reinterpolate = 40
# adaptively update the spring constant
NEBparams.adjustk_freq = 20
# change the image density
NEBparams.image_density = 5
# change the iteration density (based on the number of images)
NEBparams.iter_density = 40

# Transition State Search parameters:
tsparams = decparams.local_connect_params.tsSearchParams
# change the maximum uphill step
tsparams.max_uphill_step = .05
# change the number of steps in the transverse minimization
tsparams.nsteps_tangent1 = 4

#
# done modifying parameters
#


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
