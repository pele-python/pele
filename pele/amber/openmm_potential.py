"""
Wrapper for OpenMM Amber potential.
     
To be consistent with GMIN, units are kcal/mol and angstroms 

Requires: 
         coords.inpcrd and coords.prmtop 
"""
from __future__ import print_function


# TODO: if BasePotential is imported after simtk imports, it gives a seg fault!! 
from .simtk.openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation
from .simtk.openmm import *
from .simtk.unit import kilocalories_per_mole, kilojoules_per_mole, nanometer, angstrom, picosecond
from .simtk.openmm.app import forcefield as openmmff

from pele.potentials import BasePotential

__all__ = ["OpenMMAmberPotential"]


class OpenMMAmberPotential(BasePotential):
    """ 
    OpenMM  
    
    V(r) = Amber 
    """

    def __init__(self, prmtopFname, inpcrdFname):  # prmtopFname, inpcrdFname ):

        self.prmtop = AmberPrmtopFile(prmtopFname)
        self.inpcrd = AmberInpcrdFile(inpcrdFname)
        # number of atoms
        self.natoms = self.prmtop.topology._numAtoms

        # todo: set up ff and simulation object  
        self.system = self.prmtop.createSystem(nonbondedMethod=openmmff.NoCutoff)  # no cutoff
        self.integrator = VerletIntegrator(0.001 * picosecond)
        self.simulation = Simulation(self.prmtop.topology, self.system, self.integrator)

        # Another way of setting up potential using just pdb file ( no prmtop )
        # pdb = PDBFile('coords.pdb')
        # forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
        # system = forcefield.createSystem(pdb.topology, nonbondedMethod=ff.NoCutoff)
        # integrator = VerletIntegrator(0.001*picoseconds)
        # simulation = Simulation(pdb.topology, system, integrator)
        # simulation.context.setPositions(pdb.positions)

        # remove units 
        self.localCoords = self.inpcrd.positions / angstrom
        self.kJtokCal = kilocalories_per_mole / kilojoules_per_mole

    # '''  ------------------------------------------------------------------- '''
    def copyToLocalCoords(self, coords):
        """ copy to local coords -- deprecated  """

        # copy to local coords         
        for i in range(self.natoms):
            self.localCoords[i] = Vec3(coords[3 * i], coords[3 * i + 1], coords[3 * i + 2])

        # '''  ------------------------------------------------------------------- '''

    def getEnergy(self, coords):
        """ returns energy in kcal/mol """

        # using unit.Quantity is 5 times faster than calling copyToLocalCoords!  
        coordinates = unit.Quantity(coords.reshape(self.natoms, 3), angstrom)
        self.simulation.context.setPositions(coordinates)

        # attach units to coordinates before computing energy 
        self.simulation.context.setPositions(coordinates)
        potE = self.simulation.context.getState(getEnergy=True).getPotentialEnergy()

        # remove units from potE and then convert to kcal/mol to be consistent with GMIN   
        ee = potE / kilojoules_per_mole / self.kJtokCal
        return float(ee)

    # '''  ------------------------------------------------------------------- '''
    def getEnergyGradient(self, coords):
        """ returns energy and gradient in kcal/mol and kcal/mol/angstrom"""

        coordinates = unit.Quantity(coords.reshape(self.natoms, 3), angstrom)
        self.simulation.context.setPositions(coordinates)

        # get pot energy 
        potE = self.simulation.context.getState(getEnergy=True).getPotentialEnergy()

        E = potE / kilojoules_per_mole / self.kJtokCal
        # get forces  
        forcee = self.simulation.context.getState(getForces=True).getForces(asNumpy=True)
        # xply to -1 to convert gradient ; divide by 10 to convert to kJ/mol/angstroms 
        grad = -forcee / ( kilojoules_per_mole / nanometer ) / 10 / self.kJtokCal  # todo - 10 is hardcoded

        # remove units before returning  
        grad = np.array(grad, dtype=float)
        return float(E), grad.reshape(-1)


'''  ------------------------------------------------------------------- '''

if __name__ == "__main__":
    # create potential for the molecule in coords.pdb
    prmtopFname = '../../examples/amber/aladipep/coords.prmtop'
    inpcrdFname = '../../examples/amber/aladipep/coords.inpcrd'
    pot = OpenMMAmberPotential(prmtopFname, inpcrdFname)

    # read coordinates from pdb file 
    from .simtk.openmm.app import pdbfile as openmmpdb

    pdb = openmmpdb.PDBFile('../../examples/amber/aladipep/coords.pdb')

    coords = pdb.getPositions() / angstrom
    coords = numpy.reshape(numpy.transpose(coords), 3 * len(coords), 1)

    # compute energy and gradients       
    e = pot.getEnergy(coords)
    print('Energy (kJ/mol) = ')
    print(e)

    e, g = pot.getEnergyGradient(coords)
    gnum = pot.NumericalDerivative(coords, eps=1e-6)

    print('Energy (kJ/mol) = ')
    print(e)
    print('Analytic Gradient = ')
    print(g[60:65])
    print('Numerical Gradient = ')
    print(gnum[60:65])

    import numpy as np

    print('Num vs Analytic Gradient =')
    print(np.max(np.abs(gnum - g)), np.max(np.abs(gnum)))
    print(np.max(np.abs(gnum - g)) / np.max(np.abs(gnum)))


# $ python openmm_potential.py
# Energy (kJ/mol) = 
# -13.2272103351
# Energy (kJ/mol) = 
# -13.2272103351
# Analytic Gradient = 
# [0.5474055271317946 -1.3862268604248615 -0.34820375884655835
# 0.557332675020795 -1.3579497787292465]
# Numerical Gradient = 
# [ 0.21390981 -1.1717879  -0.09397455  0.1569264  -1.22877131]
# Num vs Analytic Gradient =
# 1.93645994676 19.2856516597
# 0.100409360334




