import ambgmin_ as GMIN
import pygmin.potentials.gminpotential as gminpot

# OpenMM 
from simtk.openmm.app import * # AmberPrmtopFile, AmberInpcrdFile, PDBFile, Simulation
from simtk.openmm import * 
from simtk.unit import picosecond
import simtk.openmm.app.forcefield as ff

#-------------------------------------------------------------------
#
# This script compares Amber energies from GMIN binding and two different ways via OpenMM  
# 
# GMIN Input files are coords.inpcrd, coords.prmtop, min.in, data (required but not used). 
#   The energy from Fortran GMIN is -21.7345926639 kcal/mol    
# 
# OpenMM calculations: 
#   1. coords.inpcrd for coordinates and coords.prmtop for ff params 
#   2. coords.pdb for coordinates and picks Amber ff params from OpenMM's built-in implementation
#  
#    Strangely the second calculation is in better agreement with GMIN energy! 
#
#-------------------------------------------------------------------

# energy from GMIN 
GMIN.initialize()   # reads coords.inpcrd and coords.prmtop 
pot = gminpot.GMINPotential(GMIN)

coords = pot.getCoords()
enerGmin = pot.getEnergy(coords)*4.184

# ----- OpenMM
 
# setup using inpcrd and prmtop 
prmtop = AmberPrmtopFile('coords.prmtop')
inpcrd = AmberInpcrdFile('coords.inpcrd')
system1 = prmtop.createSystem(nonbondedMethod=ff.NoCutoff )
integrator1 = VerletIntegrator(0.001*picosecond)
simulation1 = Simulation(prmtop.topology, system1, integrator1)
simulation1.context.setPositions(inpcrd.positions)
# get energy 
ener1 = simulation1.context.getState(getEnergy=True).getPotentialEnergy()

# setup using pdb and built-in amber ff
pdb = PDBFile('coords.pdb')
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
system2 = forcefield.createSystem(pdb.topology, nonbondedMethod=ff.NoCutoff)
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

