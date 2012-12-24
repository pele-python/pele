import ambgmin_ as GMIN
import pygmin.potentials.gminpotential as gminpot
import numpy as np
import pygmin.basinhopping as bh
from pygmin.takestep import displace
from pygmin import defaults

# OpenMM 
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.forcefield as ff
from sys import stdout

# energy from GMIN 
GMIN.initialize()   
pot = gminpot.GMINPotential(GMIN)

coords = pot.getCoords()
print pot.getEnergy(coords)*4.184

# energy from OpenMM (using inpcrd and prmtop) 
prmtop = AmberPrmtopFile('coords.prmtop')
inpcrd = AmberInpcrdFile('coords.inpcrd')
system = prmtop.createSystem(nonbondedMethod=ff.NoCutoff )
integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
#
# energy from OpenMM (pdb and in-built amber ff) 
#pdb = PDBFile('coords.pdb')
#forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
#system = forcefield.createSystem(pdb.topology, nonbondedMethod=ff.NoCutoff)
#integrator = VerletIntegrator(0.001*picoseconds)
#simulation = Simulation(pdb.topology, system, integrator)
#simulation.context.setPositions(pdb.positions)

# get OpenMM energy 
ener = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print ener

simulation.reporters.append(StateDataReporter(stdout, 1, step=True, potentialEnergy=True, temperature=True))
simulation.step(1)

print "--------------"

# BH continues 
step = displace.RandomDisplacement(stepsize=0.7)

opt = bh.BasinHopping(coords, pot, takeStep=step)

opt.quenchParameters['tol'] = 1e+10 # 1e-4
opt.run(3)

# some visualization
#try: 
#    import pygmin.utils.pymolwrapper as pym
#    pym.start()
#    pym.draw_spheres(opt.coords, "A", 1)
#except:
#    print "Could not draw using pymol, skipping this step"
#import pygmin.printing.print_atoms_xyz as pr
#
#pr.printAtomsXYZ(open("final.xyz", "w"), opt.coords)


