
from pygmin.systems.amberSystem import AMBERSystem_OpenMM

# create a new amber system
sys  = AMBERSystem_OpenMM('coords.prmtop', 'coords.inpcrd')

# -- check a single conformation 
coords = sys.get_random_configuration()
print 'Check result = ', sys.check_cistrans(coords)

# -- check minima and TS in an existing database
 
dbcurr = sys.create_database(db="aladipep.db")

print 'Checking minima .. '
ct = 0 
for minimum in dbcurr.minima():
#    print minimum.coords     
    if sys.check_cistrans(minimum.coords):
        print 'PASS', minimum._id, minimum.energy
        ct = ct + 1 
        dbcurr.removeMinimum(minimum)
    else: 
        print 'FAIL', minimum._id, minimum.energy
        
    print '------------'
        
print 'Number of minima deleted = ', ct 

print 'Checking transition states .. '
ct = 0 
print len(dbcurr.transition_states())
for ts in dbcurr.transition_states() :
    if sys.check_cistrans(ts.coords ):
            print 'PASS', ts._id, ts.energy
            ct = ct + 1 
    #        dbcurr.removeTS(ts) # not implemented yet 
    else: 
            print 'FAIL', ts._id, ts.energy            
            
    print '------------'
    
print 'Number of TS deleted = ', ct 
            










