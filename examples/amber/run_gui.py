from pygmin.systems.amberSystem import AMBERSystem

# create new amber system
print '----------------------------------'
print 'GMIN POTENTIAL' 
sys   = AMBERSystem('coords.prmtop', 'coords.inpcrd')        

#start the gui 
from pygmin.gui import run as gr    
gr.run_gui(sys)
