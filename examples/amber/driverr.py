
from pygmin.systems.amberSystem import AMBERSystem_GMIN, AMBERSystem_OpenMM 

# create new amber system
print '----------------------------------'
print 'GMIN POTENTIAL' 
sys   = AMBERSystem_GMIN('coords.prmtop', 'coords.inpcrd')        
sys.test_potential('coords.pdb')

#print 'OPENmm POTENTIAL' 
#sysOpenMM  = AMBERSystem_OpenMM('coords.prmtop', 'coords.inpcrd')
#sysOpenMM.test_potential('coords.pdb')

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
from pygmin.gui import run as gr    
gr.run_gui(sys, db="aladipep.db")

# ------ Test potential 
sys.test_potential('coords.pdb')
    
# ------ BH 
sys.test_BH(dbcurr)

# ------- Connect runs 
sys.test_connect(dbcurr)  
    
# ------- Disconn graph  
sys.test_disconn_graph(dbcurr)  
    
# ------- Test mindist  
sys.test_mindist(dbcurr)
    








