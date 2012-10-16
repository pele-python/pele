import numpy as np
from copy import copy
from pygmin.potentials.potential import potential as BasePotential

from simtk.openmm.app import * # AmberPrmtopFile, AmberInpcrdFile, PDBFile, Simulation
from simtk.openmm import * 
from simtk.unit import nanometers, picosecond, kilojoule_per_mole
import simtk.openmm.app.forcefield as ff

__all__ = ["openmmAmberPotential"]
          

class openmmAmberPotential(BasePotential):
    """
    openmmAmberPotential
    """
    def __init__(self, pdbfname ):
                
        # setup using pdb and built-in amber ff
        self.pdb = PDBFile(pdbfname)
        self.coords= self.pdb.getPositions(asNumpy=True)
        
        self.natoms = np.size(self.coords, 0)
        print self.natoms 
         
        forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
        system = forcefield.createSystem(self.pdb.topology, nonbondedMethod=ff.NoCutoff)
        integrator = VerletIntegrator(0.001*picosecond)
        self.simulation = Simulation(self.pdb.topology, system, integrator)
                                
#        self.simulation.context.setPositions(self.pdb.positions) # okay 

#        self.coords= stripUnits(self.coords) 
#        self.simulation.context.setPositions(self.coords) # works only if coords does not have units  
        
    def getEnergy(self, coordsVec ):
        
        # copy to data structure which can be passed to openmm  
        coordsLoc = [] ;         
        for i in range(coordsVec.size/3):
            coordsLoc.append(Vec3(coordsVec[3*i],coordsVec[3*i+1],coordsVec[3*i+2] ))
        
        self.simulation.context.setPositions(coordsLoc*nanometers)           

        state = self.simulation.context.getState(getEnergy=True)
        E = state.getPotentialEnergy()
        return E/ kilojoule_per_mole # to return a float 
        
#    def getEnergyGradient(self):
#        #do internal energies first
#        E = 0.
#        grad = np.zeros(self.natoms)
#        return E, grad
        
if __name__ == "__main__":
    
    pdbfname = '/home/ss2029/WORK/PyGMIN/examples/OpenMM/compare_with_ambgmin/coords.pdb'
    pot = openmmAmberPotential( pdbfname )
    
    coordsVec = np.random.uniform(-1,1,pot.natoms*3)*2
 
    # get coords from pdb file 
    pdb1 = PDBFile(pdbfname)
    coords1= pdb1.getPositions(asNumpy=True)/nanometers     
    ct=0 
    for i in coords1:
        for j in i:            
            coordsVec[ct]=j
            ct=ct+1  
            
    e = pot.getEnergy(coordsVec)    
    print "energy ", e
            
    print "numerical gradient"
#    ret = pot.getEnergyGradientNumerical(coordsVec)
    ret = pot.getEnergyGradient(coordsVec) # effectively numerical 
#    print ret[1]

########################################################
    print "trying a quench"
    
    from pygmin.optimize.quench import quench, cg , fire 

# lbfgs 
#    ret = quench( coordsVec, pot.getEnergyGradient, iprint=-1 , tol = 1e-3, nsteps=1) 
     # core dump! 

# cg  
#    ret = cg( coordsVec, pot.getEnergyGradient) 
    # runtime error -- ValueError: The truth value of an array with more than ...
    
# fire   
    ret = fire( coordsVec, pot.getEnergyGradient, tol = 1e-3, nsteps=1000) # ValueError: The truth value of an array with more than ...
    # works but after 1000 iterations gives an energy of -90.9378267921 higher than initial energy of -90.9364375726!
    
    print "energy ", ret[1]
    print "rms gradient", ret[2]
    print "number of function calls", ret[3]    
        



