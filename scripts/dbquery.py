import sys
from pygmin.storage.database import Database
import os
import getopt


def printEnergies(database):
    for min1 in database.minima():
        print "energy ", min1.energy, min1._id

def printTS(database):
    for ts in database.transition_states():
        m1 = ts.minimum1
        m2 = ts.minimum2
        dist = database.getDistance(m1, m2)
        print "ts %11.3f eigenvalue %9.3f connects %4d %4d energies %11.3f %11.3f distance %.3f" %(
        ts.energy, ts.eigenval, m1._id, m2._id, m1.energy, m2.energy, dist)
        #print "ts ", ts.energy, "eigenvalue", ts.eigenval, "connects", m1._id, m2._id, "energies", m1.energy, m2.energy, "distance", dist

def printDist(database):
    for distance in database.distances():
        m1 = distance.minimum1
        m2 = distance.minimum2
        print "distance: %4d %4d distance %.3f" %(m1._id, m2._id, distance.dist)
        #print "ts ", ts.energy, "eigenvalue", ts.eigenval, "connects", m1._id, m2._id, "energies", m1.energy, m2.energy, "distance", dist

def printCoords(id, database):
    min1 = None
    for m in database.minima():
        if m._id == id:
            min1=m
            break
    if min1 is None:
        return
    print "#energy", min1.energy
    coords = min1.coords.reshape([-1,3])
    for i in range(len(coords[:,0])):
        print "%f %f %f" % tuple(coords[i,:])

if __name__ == "__main__":
    usage = "python %s [options] dbfile\n" % sys.argv[0]
    usage += """Print minima and transition state information from the database 
        -h, --help : print this help info
        --dist     : print distances as well
        --dump-coords=id : dump the coords of minima with id=id 
    """

    printe = True
    printts = True
    printdist = False
    dumpcoords = False
    opts, args = getopt.gnu_getopt(sys.argv[1:], "h", ["help", "dist", "dump-coords="])
    for o, a in opts:
        if o == "--dist":
            printdist = True
        elif o in ("-h", "--help"):
            print usage
            exit(1)
        elif o == "--dump-coords":
            dumpcoords=True
            id = int(a)
        else:
            assert False, "unhandled option"

    if len(args) == 0:
        print usage
        exit()
    else:
        dbfile = args[0]
        
    database = Database(dbfile)
    
    if dumpcoords:
        printCoords(id, database)
        exit(1)
    
    if printe: printEnergies(database)
    if printts: printTS(database)
    if printdist: printDist(database)
    
