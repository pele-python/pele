"""
Wrapper for OpenMM Amber potential     
"""

from pygmin.potentials import BasePotential
import ambgmin_ as GMIN 
from pygmin.potentials.gminpotential import GMINPotential

# OpenMM - just read prmtop and crd file 
from simtk.openmm.app import AmberPrmtopFile, AmberInpcrdFile
from simtk.unit import angstrom as openmm_angstrom

__all__ = ["GMINAmberPotential"]

class GMINAmberPotential(GMINPotential):
    """ 
    GMIN amber potential  
    
    V(r) = Amber ff 
 
    """
    def __init__(self, prmtopFname, inpcrdFname  ): # prmtopFname, inpcrdFname ):
        # reads coords.inpcrd , coords.prmtop , min.in and data 
        #  - fnames hard coded (todo)
        super(GMINAmberPotential, self).__init__(GMIN)
#        self.potentialLocal = gminpot.GMINPotential(GMIN)
        GMIN.initialize()                
        
        self.prmtop = AmberPrmtopFile( prmtopFname )
        self.inpcrd = AmberInpcrdFile( inpcrdFname ) 
        # number of atoms
        self.natoms = self.prmtop.topology._numAtoms
                            
# '''  ------------------------------------------------------------------- '''
#    def getEnergy(self, coords):
#        enerGmin = self.potentialLocal.getEnergy(coords)
#        
#        return enerGmin  
#
##'''  ------------------------------------------------------------------- '''
#    def getEnergyGradient(self, coords):        
#        E,gminEGrad = self.potentialLocal.getEnergyGradient(coords)                                 
#        return E, gminEGrad   
    

if __name__ == "__main__":
    print 'test via amberSystem.py' 
        