Loading an Optim Database in PELE
=================================

Input files: 
coords.inpcrd coords.prmtop 
min.data  points.min points.ts  ts.data
coordsModTerm.pdb (optional, needed for connect runs) 

$ python ../../../../scripts/optim_database_converter.py 
reading from points.min
read 22 minimum coordinates of length 66
reading from min.data
--->finished loading 22 minima
reading from points.ts
reading from ts.data
--->finished loading 35 transition states

# creates optimdb.sqllite 

$ more run_gui.py 
from pele.amber.amberSystem import AMBERSystem 
from pele.gui import run as gr    
sysAmb  = AMBERSystem('coords.prmtop', 'coords.inpcrd')
gr.run_gui(sysAmb, db="optimdb.sqlite")

$ python run_gui.py   





