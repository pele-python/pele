from pele.amber.openmm_potential import OpenMMAmberPotential 
from simtk.unit import   kilocalories_per_mole, kilojoules_per_mole, nanometer, angstrom, picosecond 
import numpy 

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

