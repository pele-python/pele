from pygmin.systems.amberSystem import AMBERSystem_GMIN

# create new amber system
print '----------------------------------'
print 'GMIN POTENTIAL' 
sys   = AMBERSystem_GMIN('coords.prmtop', 'coords.inpcrd')        

#start the gui 
from pygmin.gui import run as gr    
gr.run_gui(sys)#, db="aladipep.db")
