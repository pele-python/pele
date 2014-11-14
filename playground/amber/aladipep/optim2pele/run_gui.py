
from pele.amber.amberSystem import AMBERSystem 
from pele.gui import run as gr    

sysAmb  = AMBERSystem('coords.prmtop', 'coords.inpcrd')
gr.run_gui(sysAmb, db="optimdb.sqlite")

