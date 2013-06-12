"""

This script compares Amber energies from GMIN binding and two different ways via OpenMM.  

GMIN Input files are coords.inpcrd, coords.prmtop and min.in. From Fortran code the energy is -21.7345926639 kcal/mol    

One of the OpenMM calculation uses coords.inpcrd for coordinates and coords.prmtop for ff params. 
The other OpenMM calc uses coords.pdb for coordinates and picks Amber ff params from OpenMM's own implementation.

Strangely the second calculation is in better agreement with GMIN energy! 

Amber system class in not used here. So this script would be a good starting point to understand how OpenMM and GMIN function calls work. 

"""

import ambgmin_ as GMIN
import pele.potentials.gminpotential as gminpot

# OpenMM 
from simtk.openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation
from simtk.openmm.app import pdbfile as openmmpdb
from simtk.openmm import * 
from simtk.unit import picosecond
import simtk.openmm.app.forcefield as openmmff
#from sys import stdout


# energy from GMIN 
GMIN.initialize()   # reads coords.inpcrd and coords.prmtop 
pot = gminpot.GMINPotential(GMIN)

coords = pot.getCoords()
enerGmin = pot.getEnergy(coords)*4.184

# ----- OpenMM
 
# setup using inpcrd and prmtop 
prmtop = AmberPrmtopFile('coords.prmtop')
inpcrd = AmberInpcrdFile('coords.inpcrd')
system1 = prmtop.createSystem(nonbondedMethod=openmmff.NoCutoff )
integrator1 = VerletIntegrator(0.001*picosecond)
simulation1 = Simulation(prmtop.topology, system1, integrator1)
simulation1.context.setPositions(inpcrd.positions)
# get energy 
ener1 = simulation1.context.getState(getEnergy=True).getPotentialEnergy()

# setup using pdb and built-in amber ff
pdb = openmmpdb.PDBFile('coords.pdb')
forcefield = openmmff.ForceField('amber99sb.xml', 'tip3p.xml')
system2 = forcefield.createSystem(pdb.topology, nonbondedMethod=openmmff.NoCutoff)
integrator2 = VerletIntegrator(0.001*picosecond)
simulation2 = Simulation(pdb.topology, system2, integrator2)
simulation2.context.setPositions(pdb.positions)
# get energy 
ener2 = simulation2.context.getState(getEnergy=True).getPotentialEnergy()

# print all energies 

print "Energies (kJ/mol)"
print "AMBGMIN        OpenMM inpcrd/prmtop   OpenMM pdb/amb99sb  "
print "-------------------------------------------------------- "
print enerGmin , ener1,  ener2 

