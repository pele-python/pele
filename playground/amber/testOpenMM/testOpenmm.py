"""

Test energies from OpenMM
-------------------------
This is just for testing the OpenMM installation. It does not test OpenMM's integration with PELE. 
Amber energies are computed in two different ways: 

1. coords.inpcrd and coords.prmtop 
2. coords.pdb  ( using built-in amber ff ) 

Usage: 

Tested with OpenMM 6 

$ python testOpenmm.py 
Energies (kJ/mol)
inpcrd/prmtop   pdb (+built-in amb99sb)  
-------------------------------------------------------- 
-55.3758515282 kJ/mol -55.342753321 kJ/mol

Last Updated: 
 23 Mar 2014 

"""
# OpenMM 
from simtk.openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation
from simtk.openmm.app import pdbfile as openmmpdb
from simtk.openmm import * 
from simtk.unit import picosecond
import simtk.openmm.app.forcefield as openmmff
#from sys import stdout

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
print "inpcrd/prmtop   pdb (+built-in amb99sb)  "
print "-------------------------------------------------------- "
print ener1,  ener2 

