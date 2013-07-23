from pygmin.amber import amberSystem   
import time  

# create new amber system
sys   = amberSystem.AMBERSystem('coords.prmtop', 'coords.inpcrd')        

# openmm potential is ~6x slower than gmin potential 

sys.test_potential('coords.pdb')

# load existing database 
from pygmin.storage import Database
dbcurr = sys.create_database(db="aladipep.db")

print "---------id, minener"

for minimum in dbcurr.minima():
    print minimum._id, minimum.energy    

print "---------id, m1_id, m2_id, tsener"
for ts in dbcurr.transition_states() :
    print ts._id, ts._minimum1_id, ts._minimum2_id,  ts.energy      
            
# create new database  
# dbcurr = sys.create_database(db=dbcurr)    

# connect to existing db 
#    sysOpenMM.create_database(db=dbcurr)    
        
# ------- TEST gui 
#from pygmin.gui import run as gr    
#gr.run_gui(sys, db="aladipep.db")

# ------ Test potential 
sys.test_potential('coords.pdb')
    
# ------ BH 

start = time.clock()
sys.test_BH(dbcurr,100)
elapsed = (time.clock() - start)
#print "time taken by BH = ", elapsed 

start = time.clock()
group_rotation_parameters = {}
group_rotation_parameters["NME"] = {}
group_rotation_parameters["NME"]["bond_atom_1"] = 8
group_rotation_parameters["NME"]["bond_atom_2"] = 14
group_rotation_parameters["NME"]["group_atoms"] = 15, 16, 17, 18, 19, 20, 21
group_rotation_parameters["NME"]["max_angle_magnitude"] = 0.3
group_rotation_parameters["NME"]["selection_probability"] = 1.0
group_rotation_parameters["ACE"] = {}
group_rotation_parameters["ACE"]["bond_atom_1"] = 8
group_rotation_parameters["ACE"]["bond_atom_2"] = 6
group_rotation_parameters["ACE"]["group_atoms"] = 0, 1, 2, 3, 4, 5, 7 
group_rotation_parameters["ACE"]["max_angle_magnitude"] = 0.3
group_rotation_parameters["ACE"]["selection_probability"] = 1.0
sys.test_BH_group_rotation(dbcurr,100, group_rotation_parameters)
elapsed2 = (time.clock() - start)
print "time taken by BH = ", elapsed
print "time taken by BH with group rotation = ", elapsed2


exit() 
# ------- Connect runs 
sys.test_connect(dbcurr)  
    
# ------- Disconn graph  
sys.test_disconn_graph(dbcurr)  
    
# ------- Test mindist  
sys.test_mindist(dbcurr)
    
