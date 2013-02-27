from pygmin.amber import amberSystem 

# create new amber system
sys   = amberSystem.AMBERSystem('coords.prmtop', 'coords.inpcrd')        

#start the gui 
from pygmin.gui import run as gr    
gr.run_gui(sys)
