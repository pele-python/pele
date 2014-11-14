from pele.amber import amberSystem
import playground.group_rotation.group_rotation as group_rotation
import pele.amber.read_amber as read_amber
import time  

# create new amber system
# It turns out that these are the only names that we can use, since we rely on 
# countatoms.f90 from AMBGMIN to determine the number of atoms for any fortran
# energy and gradient calls.
system   = amberSystem.AMBERSystem('coords.prmtop', 'coords.inpcrd')                

# openmm potential is ~6x slower than gmin potential 

#system.test_potential('coords.pdb')

# load existing database 
from pele.storage import Database
#dbcurr = system.create_database(db="test.db")
#dbcurr2 = system.create_database(db="test2.db")
dbcurr = system.create_database()
dbcurr2 = system.create_database()
print "---------id, minener"

for minimum in dbcurr.minima():
    print minimum._id, minimum.energy    

print "---------id, m1_id, m2_id, tsener"
for ts in dbcurr.transition_states() :
    print ts._id, ts._minimum1_id, ts._minimum2_id,  ts.energy      
            
# create new database  
# dbcurr = system.create_database(db=dbcurr)    

# connect to existing db 
#    sysOpenMM.create_database(db=dbcurr)    
        
# ------- TEST gui 
#from pele.gui import run as gr    
#gr.run_gui(sys, db="aladipep.db")

# ------ Test potential 
#system.test_potential('coords.pdb')
    
# ------ BH 

time.clock()
system.test_BH(dbcurr,100)
elapsed = time.clock()
#print "time taken by BH = ", elapsed 

time.clock()
group_rotation_parameters = read_amber.default_parameters(system.prmtop_name)
system.test_BH_group_rotation(dbcurr2,100, group_rotation_parameters)
elapsed2 = time.clock()
print "time taken by BH = ", elapsed
print "time taken by BH with group rotation = ", elapsed2
exit() 

# ------- Connect runs 
system.test_connect(dbcurr)  
    
# ------- Disconn graph  
system.test_disconn_graph(dbcurr)  
    
# ------- Test mindist  
system.test_mindist(dbcurr)
    
