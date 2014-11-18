from pele.storage.database import Database, Minimum, TransitionState
import time
import numpy as np
from simtk.openmm.app import AmberPrmtopFile

""" creates a sqlite database from min.data, ts.data, extractedmin, extractedts 

Input: 
    coords.prmtop ( for number of atoms ) 
    min.data, ts.data, extractedmin, extractedts    
    
Output: 
    storage.sqlite 

"""

# determine number of atoms from prmtop 
prmtop = AmberPrmtopFile( 'coords.prmtop' )  # TOSET 
natoms = prmtop.topology._numAtoms

# open database 
db = Database(db="storage.sqlite")   # TOSET 

def read_coords(filee):
    coords = np.zeros(3*natoms)
    for i in xrange(natoms):
        x = filee.readline().split()
        coords[i*3:i*3+3] = [float(y) for y in x]
    return coords 

# counter to keep track of added minima
mini=1
minima={}

# for timing
tt = t0 = time.time()

# open coordinate file
fcoords = open("extractedmin")   # TOSET 

print "Reading minima"
# loop over all lines in min
for line in open("min.data"):     # TOSET 
    coords = read_coords(fcoords) 
    energy, frequency, pgorder, itx, ity, itz = line.split()    
    min1 = Minimum(float(energy), coords)    
    db.session.add(min1)
    minima[mini]=min1
    if(mini%10000 == 0):
        print "commiting the next 10000 minima, %d in total"%(mini)
        db.session.commit()
    mini+=1
           

print "%.1f seconds"%(time.time() - tt)
tt = time.time()
print "Commiting changes to database"
db.session.commit()
print "%.1f seconds"%(time.time() - tt)
tt = time.time()


print "Reading transition states"
fcoords = open("extractedts")    # TOSET 
tsi=1 
for line in open("ts.data"):      # TOSET 
    coords = read_coords(fcoords) 
    energy, frequency, pgorder, min1, min2, itx, ity, itz = line.split()
    ts = TransitionState(float(energy), coords, minima[int(min1)], minima[int(min2)])
    db.session.add(ts)
    #db.addTransitionState(float(energy), None, minima[int(min1)], minima[int(min2)], commit=False)
    if(tsi%10000 == 0):
        print "commiting the next 10000 transition states, %d in total"%(tsi)
        db.session.commit()
    tsi+=1
print "%.1f seconds"%(time.time() - tt)
tt = time.time()
print "Commiting changes to database"
db.session.commit()
print "%.1f seconds"%(time.time() - tt)
print "Done after %.1f seconds"%(time.time() - t0)


