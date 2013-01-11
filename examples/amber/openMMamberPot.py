"""
Wrapper for OpenMM Amber potential     
  Requires: 
         coords.inpcrd and coords.prmtop 
"""

# TODO: sort out units. to be consistent with GMIN make it kcal/mol and angstroms 

# TODO: if BasePotential is imported after simtk imports, it gives a seg fault!! 
from pygmin.potentials import BasePotential

from simtk.openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation
from simtk.openmm import * 
from simtk.unit import  * 
import simtk.openmm.app.forcefield as openmmff

class amberPot(BasePotential):
    """ 
    OpenMM  
    
    V(r) = Amber 
    """
    def __init__(self, prmtopFname, inpcrdFname  ): # prmtopFname, inpcrdFname ):

        self.prmtop = AmberPrmtopFile( prmtopFname )
        self.inpcrd = AmberInpcrdFile( inpcrdFname ) 
        # number of atoms
        self.natoms = self.prmtop.topology._numAtoms
        
        # todo: set up ff and simulation object  
        self.system = self.prmtop.createSystem(nonbondedMethod=openmmff.NoCutoff ) # no cutoff
        self.integrator = VerletIntegrator(0.001*picosecond)
        self.simulation = Simulation(self.prmtop.topology, self.system, self.integrator)
        
        # remove units 
        self.localCoords = self.inpcrd.positions/angstrom 
            
# '''  ------------------------------------------------------------------- '''
    def copyToLocalCoords(self, coords):
        # copy to local coords         
        for i in range(self.natoms):
            self.localCoords[i] = Vec3(coords[3*i], coords[3*i+1] , coords[3*i+2] )              

# '''  ------------------------------------------------------------------- '''
    def getEnergy(self, coords):
        # copy to local coords 
        self.copyToLocalCoords(coords)
                
        # attach units to coordinates before computing energy 
        self.simulation.context.setPositions(self.localCoords * angstrom)
        potE = self.simulation.context.getState(getEnergy=True).getPotentialEnergy()

        # remove units before returning   todo 
        return potE / kilojoules_per_mole 

#'''  ------------------------------------------------------------------- '''
    def getEnergyGradient(self, coords):
        # copy to local coords 
        self.copyToLocalCoords(coords)
        
        # attach units to coordinates  
        self.simulation.context.setPositions(self.localCoords * angstrom)
        # get pot energy 
        potE = self.simulation.context.getState(getEnergy=True).getPotentialEnergy()
        E = potE / kilojoules_per_mole 
        # get forces  
        forcee = self.simulation.context.getState(getForces=True).getForces(asNumpy=True)
        # xply to -1 to convert gradient ; divide by 10 to convert to kJ/mol/angstroms 
        grad = -forcee / ( kilojoules_per_mole / nanometer ) / 10 # todo - 10 is hardcoded   
                
        # remove units before returning   
        return E, grad.reshape(-1)  
'''  ------------------------------------------------------------------- '''

if __name__ == "__main__":
    
    # read one conformation from pdb file 
    from simtk.openmm.app import pdbfile as openmmpdb
    pdb = openmmpdb.PDBFile('coords.pdb')
    
    coords = pdb.getPositions() / angstrom   
    coords = numpy.reshape(numpy.transpose(coords), 3*len(coords), 1)

    pot = amberPot()    
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
    
    
