import ambgmin_ as GMIN
import pele.potentials.gminpotential as gminpot

# OpenMM 
from simtk.openmm.app import * # AmberPrmtopFile, AmberInpcrdFile, PDBFile, Simulation
from simtk.openmm import * 
from simtk.unit import * 
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
E,gminEGrad = pot.getEnergyGradient(coords) 

# ----- OpenMM
 
# setup using inpcrd and prmtop 
prmtop = AmberPrmtopFile('coords.prmtop')
inpcrd = AmberInpcrdFile('coords.inpcrd')
system1 = prmtop.createSystem(nonbondedMethod=ff.NoCutoff )
integrator1 = VerletIntegrator(0.001*picoseconds)
simulation1 = Simulation(prmtop.topology, system1, integrator1)
simulation1.context.setPositions(inpcrd.positions)
# get energy and forces
state1 = simulation1.context.getState(getEnergy=True)
ener1 = state1.getPotentialEnergy()
state1 = simulation1.context.getState(getForces=True)
force1 = state1.getForces(asNumpy=True) 

# setup using pdb and built-in amber ff
pdb = PDBFile('coords.pdb')
forcefield = ForceField('amber99sb.xml') # , 'tip3p.xml')
system2 = forcefield.createSystem(pdb.topology, nonbondedMethod=ff.NoCutoff)
#system2 = forcefield.createSystem(prmtop.topology, nonbondedMethod=ff.NoCutoff)
integrator2 = VerletIntegrator(0.001*picoseconds)
simulation2 = Simulation(pdb.topology, system2, integrator2)

#simulation2.context.setPositions(pdb.positions)
simulation2.context.setPositions(inpcrd.positions)
# get energy and forces
state2 = simulation2.context.getState(getEnergy=True)
ener2 = state2.getPotentialEnergy()
state2 = simulation2.context.getState(getForces=True)
force2 = state2.getForces(asNumpy=True) 

# print all energies 
print "Energies (kJ/mol)"
print "AMBGMIN        OpenMM inpcrd/prmtop   OpenMM pdb/amb99sb  "
print "-------------------------------------------------------- "
print enerGmin , ener1,  ener2 

print "\nForces"
print force1[0,:]
print force2[0,:]
print gminEGrad[0:3]*4.184*10

print "\ncoordinates"
print inpcrd.positions[0] 
print pdb.positions[0] 

#Output: 
#Energies (kJ/mol)
#OpenMM inpcrd/prmtop   OpenMM pdb/amb99sb  
#-------------------------------------------------------- 
#-90.9428083928 kJ/mol -90.9364375726 kJ/mol
#
#Forces
#[ 0.00164192 -0.00081001  0.00142779] kJ/(nm mol)
#[-3.93361139 -0.87415225  2.41254514] kJ/(nm mol)
#
#coordinates
#(-0.8823998, 2.6018494, -1.8615158) A
#(-0.0882, 0.2602, -0.18620000000000003) nm

