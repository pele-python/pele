from pygmin.storage.database import Database, Minimum, TransitionState
import time
import numpy as np
from pygmin.potentials.lj import LJ

natoms = 38
pot = LJ()

# open database
db = Database(db="storagetmp.sqlite")

# counter to keep track of added minima
mini=1
minima={}

# for timing
tt = t0 = time.time()

# open coordinate file
fcoords = open("extractedmin")

print "Reading minima"
# loop over all lines in min
for line in open("min.data"):
    coords = np.zeros(3*natoms)
    for i in xrange(natoms):
        x = fcoords.readline().split()
        coords[i*3:i*3+3] = [float(y) for y in x] 
    energy, frequency, pgorder, itx, ity, itz = line.split()    
    print energy, pot.getEnergy(coords)
    min1 = Minimum(float(energy), coords)    
    db.session.add(min1)
    minima[mini]=min1
    mini+=1    

print "%.1f seconds"%(time.time() - tt)
tt = time.time()
print "Commiting changes to database"
db.session.commit()
print "%.1f seconds"%(time.time() - tt)
tt = time.time()


print "Reading transition states"
for line in open("ts.data"):
    energy, frequency, pgorder, min1, min2, itx, ity, itz = line.split()
    ts = TransitionState(float(energy), None, minima[int(min1)], minima[int(min2)])
    db.session.add(ts)
    #db.addTransitionState(float(energy), None, minima[int(min1)], minima[int(min2)], commit=False)
print "%.1f seconds"%(time.time() - tt)
tt = time.time()
print "Commiting changes to database"
db.session.commit()
print "%.1f seconds"%(time.time() - tt)
print "Done after %.1f seconds"%(time.time() - t0)

from  pygmin.NEB.graph import Graph
print "generating graph"
tt = time.time()
graph = Graph(db)
print "%.1f seconds"%(time.time() - tt)

#print "printing minima"
#for i in db.minima():
#    print i