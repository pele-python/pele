"""
example for how to use double ended connect to connect the minima in an existing database

we will use as an example system the Lennard-Jones cluster with a small number of atoms.
Since we don't already have a database, for this example we'll build a small one using
basinhopping"
"""
import numpy as np
import logging

from pele.systems import LJCluster
from pele.utils.disconnectivity_graph import DisconnectivityGraph, database2graph


natoms = 13
system = LJCluster(natoms)

use_existing_database = False
if use_existing_database:
    db = system.create_database("lj13.sqlite", createdb=False)
else:
    # build a small database using basinhopping
    print "building a small database using basinhopping"""
    db = system.create_database()
    bh = system.get_basinhopping(database=db, outstream=None)
    bh.run(20)

print "starting with a database of", len(db.minima()), "minima"

# turn of status printing for the connect run
# first use the logging module to turn off the status messages 
logger = logging.getLogger("pele.connect")
logger.setLevel("WARNING")


# connect all minima to the lowest minimum
print "now connecting all the minima to the lowest energy minimum"
m1 = db.minima()[0]
for m2 in db.minima()[1:]:
    print "    connecting minima with id's", m1._id, m2._id
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
