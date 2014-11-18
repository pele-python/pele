from simtk.unit import   kilocalories_per_mole, kilojoules_per_mole, nanometer, angstrom, picosecond 
import numpy 
# --- Test  OpenMMAmberPotential 
from pele.amber.openmm_potential import OpenMMAmberPotential 

# create potential for the molecule in coords.pdb
prmtopFname = '../aladipep/coords.prmtop' 
inpcrdFname	= '../aladipep/coords.inpcrd' 
pot = OpenMMAmberPotential(prmtopFname, inpcrdFname)  

# read coordinates from pdb file 
from simtk.openmm.app import pdbfile as openmmpdb
pdb = openmmpdb.PDBFile('../aladipep/coords.pdb')
    
coords = pdb.getPositions() / angstrom   
coords = numpy.reshape(numpy.transpose(coords), 3*len(coords), 1)

# compute energy and gradients   	
e = pot.getEnergy(coords)
print 'Energy (kJ/mol) = '
print e
    
e, g = pot.getEnergyGradient(coords)
gnum = pot.NumericalDerivative(coords, eps=1e-6)

print 'Energy (kJ/mol) = '
print e 
print 'Analytic Gradient = '
print g[60:65] 
print 'Numerical Gradient = '
print gnum[60:65] 

import numpy as np 
print 'Num vs Analytic Gradient =' 
print np.max(np.abs(gnum-g)), np.max(np.abs(gnum))
print np.max(np.abs(gnum-g)) / np.max(np.abs(gnum))

# --- Test  AMBERSystem class 
from pele.amber.amberSystem import AMBERSystem 

# create new amber system    
sysAmb  = AMBERSystem('../aladipep/coords.prmtop', '../aladipep/coords.inpcrd')
    
# load existing database 
from pele.storage import Database
dbcurr = Database(db="../aladipep/aladipep.db")
                            
# ------ Test potential 
print 'testing potential in ambSystem' 
sysAmb.test_potential("../aladipep/coords.pdb")
    
# ------ BH 
print 'testing BH' 
nsteps = 1
sysAmb.test_BH(dbcurr, nsteps)

for minimum in dbcurr.minima():
    print minimum._id, minimum.energy    

# -- test connect 
dbcurr = Database(db="aladipep.db")
sysAmb.test_connect(dbcurr) 


#print "---------id, m1_id, m2_id, tsener"
#for ts in dbcurr.transition_states() :
#    print ts._id, ts._minimum1_id, ts._minimum2_id,  ts.energy      
            
# connect to existing db 
#    sysOpenMM.create_database(db=dbcurr)   



# ------- TEST gui 
#from pele.gui import run as gr    
#gr.run_gui(sysAmb, db="aladipep.db")











