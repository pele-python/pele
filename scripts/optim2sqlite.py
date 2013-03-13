"""
Creates a database from 
  min.data 
  extractedmin 
  ts.data 
  extractedts 
  
  Usage: 
    run in folder containing above files 
"""

from pygmin.storage.database import Database, Minimum, TransitionState
import time
import numpy as np
import os 

def read_coords(file):
    coords = np.zeros(3*natoms)
    for i in xrange(natoms):
        x = file.readline().split()
        coords[i*3:i*3+3] = [float(y) for y in x]
    return coords 

# ---- main 

natoms = input('Enter number of atoms: ') 

# open database
print 'database file name = storage.sqlite' 

if os.path.exists('storage.sqlite'):
    db = Database(db="storage.sqlite")
    print 'database file exists'
    if hasattr(db,'minima'):
        print ' number of minima =  ', len(db.minima())
        
    if hasattr(db,'ts'):
        print  ' number of ts = ', len(db.ts() ) 

    print '\n Note that minima and ts will be appended to existing set without testing for duplicates \n'
else:
    db = Database(db="storage.sqlite")

# --- Minima 
# counter to keep track of added minima
ct=1
minima={}

# for timing
tt = t0 = time.time()

if os.path.exists('extractedmin') and os.path.exists('min.data'):
    print 'Loading minima from extractedmin and min.data ..'
    # open coordinate file
    fcoords = open("extractedmin")
    
    print "Reading minima"
    # loop over all lines in min
    for line in open("min.data"):
        coords = read_coords(fcoords) 
        energy, frequency, pgorder, itx, ity, itz = line.split()    
        min1 = Minimum(float(energy), coords)    
        min1.fvib = frequency 
        min1.pgorder = pgorder  
    
        db.session.add(min1)
        minima[ct]=min1
        if(ct%10000 == 0):
            print "commiting the next 10000 minima, %d in total"%(ct)
            db.session.commit()
        ct = ct + 1 
               
    
    print "Number of minima added = ", ct-1 
    print "%.1f seconds"%(time.time() - tt)
    tt = time.time()
    print "Commiting changes to database"
    db.session.commit()
    print "%.1f seconds"%(time.time() - tt)
    tt = time.time()
else:
    print 'Minima not added. extractedmin/min.data not found.'
    
# --- TS  

if os.path.exists('extractedts') and os.path.exists('ts.data'):

    print '\nLoading ts from extractedts and ts.data ..'
    
    fcoords = open("extractedts")
    tsi=1
    for line in open("ts.data"):
        coords = read_coords(fcoords) 
        energy, frequency, pgorder, min1, min2, itx, ity, itz = line.split()
        ts = TransitionState(float(energy), coords, minima[int(min1)], minima[int(min2)])
        ts.fvib = frequency 
        ts.pgorder = pgorder
    
        db.session.add(ts)
        #db.addTransitionState(float(energy), None, minima[int(min1)], minima[int(min2)], commit=False)
        if(tsi%10000 == 0):
            print "commiting the next 10000 transition states, %d in total"%(tsi)
            db.session.commit()
        tsi+=1

    print "Number of ts added = ", tsi-1  
    print "%.1f seconds"%(time.time() - tt)
    tt = time.time()
    print "Commiting changes to database"
    db.session.commit()
    print "%.1f seconds"%(time.time() - tt)
    print "Done after %.1f seconds"%(time.time() - t0)
else:
    print '\nTS not added. extractedts/ts.data not found. '
    