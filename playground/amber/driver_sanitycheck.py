from pele.amber import amberSystem as amb

# create a new amber system and load database to be pruned 
sys    = amb.AMBERSystem('coords.prmtop', 'coords.inpcrd')
dbcurr = sys.create_database(db="aladipep.db")

print 'Collecting minima to delete .. '
listTODel = [] 

for minimum in dbcurr.minima():
    testOutCome1 = sys.check_cistrans(minimum.coords) 
    testOutCome2 = sys.check_CAchirality(minimum.coords) 
    if testOutCome1 and testOutCome2:
        print 'PASS', minimum._id, minimum.energy
    else: 
        listTODel.append(minimum)
        print 'FAIL', minimum._id, minimum.energy
    print '------------'
        
print 'Number of minima to be deleted = ', len(listTODel) 

# now delete 
for minn in listTODel:
    dbcurr.removeMinimum(minn)


#print 'Checking transition states .. '
#ct = 0 
#print len(dbcurr.transition_states())
#for ts in dbcurr.transition_states() :
#    if sys.check_cistrans(ts.coords ):
#            print 'PASS', ts._id, ts.energy
#            ct = ct + 1 
#    #        dbcurr.removeTS(ts) # not implemented yet 
#    else: 
#            print 'FAIL', ts._id, ts.energy            
#            
#    print '------------'
#    
#print 'Number of TS deleted = ', ct 
