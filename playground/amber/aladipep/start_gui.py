from pele.amber import amberSystem 

# create new amber system
system   = amberSystem.AMBERSystem('coords.prmtop', 'coords.inpcrd')
database = system.create_database("optimdb.sqlite")    

#start the gui 
from pele.gui import run as gr    
gr.run_gui(system, db=database)
