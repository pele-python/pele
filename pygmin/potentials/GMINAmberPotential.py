"""
Wrapper for OpenMM Amber potential     
"""

from pygmin.potentials import BasePotential
import ambgmin_ as GMIN 
import pygmin.potentials.gminpotential as gminpot

# OpenMM - just read prmtop and crd file 
from simtk.openmm.app import AmberPrmtopFile, AmberInpcrdFile
from simtk.unit import angstrom as openmm_angstrom
from simtk.openmm import Vec3 
from simtk.unit import   kilocalories_per_mole, kilojoules_per_mole, nanometer, angstrom, picosecond 

__all__ = ["GMINAmberPotential"]

class GMINAmberPotential(BasePotential):
    """ 
    GMIN amber potential  
    
    V(r) = Amber ff 
 
    """
    def __init__(self, prmtopFname, inpcrdFname  ): # prmtopFname, inpcrdFname ):
        # reads coords.inpcrd , coords.prmtop , min.in and data 
        #  - fnames hard coded (todo)
        GMIN.initialize()                          
        self.potentialLocal = gminpot.GMINPotential(GMIN)
        
        self.prmtop = AmberPrmtopFile( prmtopFname )
        self.inpcrd = AmberInpcrdFile( inpcrdFname ) 
        # number of atoms
        self.natoms = self.prmtop.topology._numAtoms
        self.localCoords = self.inpcrd.positions/angstrom

# '''  ------------------------------------------------------------------- '''                            
    def copyToLocalCoords(self, coords):
        """ copy to local coords """
                
        # copy to local coords         
        for i in range(self.natoms):
            self.localCoords[i] = Vec3(coords[3*i], coords[3*i+1] , coords[3*i+2] )              
                            
# '''  ------------------------------------------------------------------- '''
    def getEnergy(self, coords):
        enerGmin = self.potentialLocal.getEnergy(coords)
        
        return enerGmin  

#'''  ------------------------------------------------------------------- '''
    def getEnergyGradient(self, coords):        
        E,gminEGrad = self.potentialLocal.getEnergyGradient(coords)                                 
        return E, gminEGrad   
    

if __name__ == "__main__":
    print 'test via amberSystem.py' 
        