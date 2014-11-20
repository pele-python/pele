"""
examples of how to do various things using the system class.  I will use
the Lennard-Jones system as an example
"""
import logging

from pele.systems import LJCluster
from pele.utils.disconnectivity_graph import DisconnectivityGraph, database2graph


# initialize the system class
natoms = 13
system = LJCluster(natoms)

# make a random configuration
coords = system.get_random_configuration()

# compute the energy of that configuration
potential = system.get_potential()
energy = potential.getEnergy(coords)
print "the energy of the random configuration is", energy

# minimize that configuration
quencher = system.get_minimizer()
ret = quencher(coords)
newcoords = ret.coords
newenergy = ret.energy
print "after quenching, the energy is", newenergy

# create a database to store the minimum in
db = system.create_database()
# note: this creates a database in memory, if you want to save the results
# you would use
# db = system.create_database("lj.sqlite")
minimum1 = db.addMinimum(newenergy, newcoords)

# get a second random minimized configuration and add it to the database
ret = system.get_random_minimized_configuration()
print "a second minimum has energy", ret[1]
e2, coords2 = ret[1], ret[0]
minimum2 = db.addMinimum(e2, coords2)

# do a basinhopping run to find the global minimum and build up the database of minima
bh = system.get_basinhopping(database=db, outstream=None)
niter = 20
bh.run(niter)
print "the lowest energy found after", niter, " basinhopping steps is", db.minima()[0].energy

# print the energies of all the minima found
print "the minima in the database have energies"
for minimum in db.minima():
    print "  ", minimum.energy

# find the minimum distance (a.k.a. mindist) between the two lowest minima
m1, m2 = db.minima()[:2]
mindist = system.get_mindist()
dist, coords1, coords2 = mindist(m1.coords, m2.coords)
print "the minimum distance between the two lowest minima is", dist

# find a connected series of minima and transition states between m1 and m2
# 
# first use the logging module to turn off the status messages 
logger = logging.getLogger("pele.connect")
logger.setLevel("WARNING")
# now create the double ended connect object
connect = system.get_double_ended_connect(m1, m2, db)
connect.connect()
mints, S, energies = connect.returnPath()
nts = (len(mints) - 1) / 2
print "\nprint found a connection with", nts, "transition states"

# connect all minima to the lowest minimum
print "now connecting all the minima to the lowest energy minimum"
from pele.landscape import ConnectManager

manager = ConnectManager(db, strategy="gmin")
for i in xrange(db.number_of_minima() - 1):
    print "connecting minima with id's", m1._id, m2._id
    m1, m2 = manager.get_connect_job()
    connect = system.get_double_ended_connect(m1, m2, db)
    connect.connect()

# print some data about the database
print "database summary:"
print "    ", len(db.minima()), "minima"
print "    ", len(db.transition_states()), "transition states"

# finally, create a disconnectivity graph from the database
print "computing and showing disconnectivity graph"
import pylab as pl

graph = database2graph(db)
dg = DisconnectivityGraph(graph, nlevels=6)
dg.calculate()
dg.plot()
pl.show()


    


