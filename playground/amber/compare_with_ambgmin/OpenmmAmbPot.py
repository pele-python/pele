import ambgmin_ as GMIN
import pele.potentials.gminpotential as gminpot
from pele.optimize import fire 

import numpy as np
from copy import copy
from pele.potentials.potential import potential as BasePotential

from simtk.openmm.app import * 
from simtk.openmm import * 
from simtk.unit import * 
import simtk.openmm.app.forcefield as ff

__all__ = ["OpenmmAmbPot"]
          
class OpenmmAmbPot(BasePotential):
    """
    Wrapper for OpenMM Amber Potential
    """
    def __init__(self, pdbfname ):
                
        # setup using pdb and built-in amber ff
        self.pdb = PDBFile(pdbfname)
        self.coords= self.pdb.getPositions(asNumpy=True)
        
        self.natoms = np.size(self.coords, 0)
        print self.natoms 

        # todo - take system setup out of constructor
        prmtop = AmberPrmtopFile('coords.prmtop')
        inpcrd = AmberInpcrdFile('coords.inpcrd')
        system = prmtop.createSystem(nonbondedMethod=ff.NoCutoff )
        integrator = VerletIntegrator(0.001*picoseconds)
        self.simulation = Simulation(prmtop.topology, system, integrator)

        # ---------------------------------------------------------------------------
                                                        
    def getEnergy(self, coordsVec ):
        
        # copy to data structure which can be passed to openmm  
        coordsLoc = [] ;         
        for i in range(coordsVec.size/3):
            coordsLoc.append(Vec3(coordsVec[3*i],coordsVec[3*i+1],coordsVec[3*i+2] ))
        
        # todo - units are hard coded  
        self.simulation.context.setPositions(coordsLoc*nanometers)           

        state = self.simulation.context.getState(getEnergy=True)
        E = state.getPotentialEnergy()
        return E/ kilojoule_per_mole # to return a float 
        # ---------------------------------------------------------------------------        

    def getEnergyGradient(self,coordsVec):        
        # copy to data structure which can be passed to openmm  
        coordsLoc = [] ;         
        for i in range(coordsVec.size/3):
            coordsLoc.append(Vec3(coordsVec[3*i],coordsVec[3*i+1],coordsVec[3*i+2] ))
        
        # todo - units are hard coded  
        self.simulation.context.setPositions(coordsLoc*nanometer)           
        state = self.simulation.context.getState(getEnergy=True, getForces=True)
        
        # remove units before returning         
        E    = state.getPotentialEnergy() / kilojoule_per_mole
        grad = state.getForces(asNumpy=True) / (kilojoule_per_mole / nanometer)
                 
        g = np.zeros(3*self.natoms)
        self.copyList2vector(g, grad ) # reshape to a 1-D array

        return E  , -g # todo: g or -g ? 
        # ---------------------------------------------------------------------------        

    def copyList2vector(self,outvec,inlist):        
        ct= 0;         
        for i in inlist:
            for j in i:
                outvec[ct] = j
                ct += 1
        # ---------------------------------------------------------------------------        

        
if __name__ == "__main__":
    
    pdbfname = 'coords.pdb'
    pot = OpenmmAmbPot( pdbfname )
    
    # get coords from pdb file 
    pdb    = PDBFile(pdbfname)   
    coords = pdb.getPositions(asNumpy=True)/nanometer # PDBFile converts coords to nm
    
    # copy coords to coordsVec      
    coordsVec = np.zeros( 3*pot.natoms, np.float64 )
    pot.copyList2vector(coordsVec, coords ) 
            
    # test getEnergy and gradients 
    e = pot.getEnergy(coordsVec)    
    print "---energy "
    print  e
            
    print "---numerical gradient"
    enum, gnum = pot.getEnergyGradientNumerical(coordsVec)
    print enum   
    print gnum[0:5]
    
    print "---openmm gradient"
    eo,go = pot.getEnergyGradient(coordsVec)  
    print eo  # energy 
    print go[0:5]    

    print "\n-----------------"    
    print "quench\n"
    
    # lbfgs 
    #  ret = quench( coordsVec, pot.getEnergyGradient, iprint=-1 , tol = 1e-3, nsteps=100) 
    # core dump! 

    # cg  
    #    ret = cg( coordsVec, pot.getEnergyGradient) 
    # runtime error -- ValueError: The truth value of an array with more than ...
    
    # fire   
    retOpmm = fire( coordsVec, pot, tol = 1e-3, nsteps=1000) 
    # works but quenched energy is higher! 
    
    print "quenched energy ", retOpmm.energy
    print "rms gradient", retOpmm.rms
    print "number of function calls", retOpmm.nfev
        
    # -------- GMIN 
    print "\n\nCompare with GMIN"    
    GMIN.initialize()   # reads coords.inpcrd and coords.prmtop 
    pot = gminpot.GMINPotential(GMIN)

    coords = pot.getCoords()
    enerGmin = pot.getEnergy(coords)
    egmin,gminEGrad = pot.getEnergyGradient(coords)
    
    retGmin = fire( coords, pot, tol = 1e-3, nsteps=1000)

    print " -- pre-quench --" 
    print "E       gmin : ", egmin 
    print "       openmm: ", eo/4.184
     
    print "grad    gmin : ", gminEGrad[0:3]
    print "       openmm: ", go[0:3]/41.84

    print " -- post-quench --" 
    print "E       gmin : ", retGmin.energy
    print "       openmm: ", retOpmm.energy/4.184     

    print "end"

