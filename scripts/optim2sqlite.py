"""
Creates a database from 
  min.data 
  extractedmin 
  ts.data 
  extractedts 
"""

from pygmin.storage.database import Database, Minimum, TransitionState
import time
import numpy as np

natoms = input('Enter number of atoms: ') 
# pot = LJ()
# open database
db = Database(db="storage.sqlite")

def read_coords(file):
    coords = np.zeros(3*natoms)
    for i in xrange(natoms):
        x = file.readline().split()
        coords[i*3:i*3+3] = [float(y) for y in x]
    return coords 


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
fcoords = open("extractedts")
tsi=1
for line in open("ts.data"):
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
