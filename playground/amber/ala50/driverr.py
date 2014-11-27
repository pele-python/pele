#from pele.amber.amberSystem import AMBERSystem_GMIN, AMBERSystem_OpenMM
from pele.amber import amberSystem
import time  

# create new amber system
#print '----------------------------------'
#print 'GMIN POTENTIAL' 
#sys   = AMBERSystem_GMIN('coords.prmtop', 'coords.inpcrd')        
#sys.test_potential('coords.pdb')

# openmm potential is ~6x slower than gmin potential 
print 'OPENmm POTENTIAL' 
sys  = AMBERSystem_OpenMM('coords.prmtop', 'coords.inpcrd')
sys.test_potential('coords.pdb')

# create new database  
from pele.storage import Database
dbcurr = sys.create_database()
                    
# ------- TEST gui 
from pele.gui import run as gr    
gr.run_gui(sys, db=dbcurr)

# ------ Test potential 
sys.test_potential('coords.pdb')
    
# ------ BH 
start = time.clock()
sys.test_BH(dbcurr)
elapsed = (time.clock() - start)
print "time taken by BH = ", elapsed 

exit() 
# ------- Connect runs 
sys.test_connect(dbcurr)  
    
# ------- Disconn graph  
sys.test_disconn_graph(dbcurr)  
    
# ------- Test mindist  
sys.test_mindist(dbcurr)
    








